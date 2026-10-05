#!/usr/bin/env python

import argparse

import numpy as np
import xarray as xr

'''
finds the 200m (by default) isobath (optionally different target depth)
given the bathymetry computed from the meshmask.
the isobath is defined on cell edges, so on the U and V grids.
saves it to a netcdf file, that can be plotted with plot_isobath.nc

Can also be imported as a library: compute_isobath(meshmask_file, depth)
runs the full extraction pipeline and returns the contours plus the mesh
arrays needed to turn them into grid indices (see create_fluxes.py)
or to plot them (see plot_isobath.py).

Usage:
source /g100/home/userexternal/ggalli00/PY_ENVS/bitsea_gg/sequence3.sh
maskfile='/g100_work/OGS_test2528/ggalli00/EIS_2026_11/qDEG_setup/MASK/meshmask_025_z125.nc'
outfile='isobath_200m.py'
depth=200.0
python define_isobath.py -m $maskfile -d $depth -o $outfile
'''


def read_meshmask(filename):
    """
    Read the NEMO meshmask.
    Bathymetry is calculated as the sum of e3t_0 over wet T cells.
    /g100_work/OGS_test2528/ggalli00/EIS_2026_11/qDEG_setup/MASK/meshmask_025_z125.nc
    """
    ds = xr.open_dataset(filename)
    tmask = ds["tmask"].isel(time=0).values.astype(bool)
    e3t = ds["e3t_0"].isel(time=0).values
    # C-grid coordinates
    glamt = ds["glamt"].isel(time=0, z_a=0).values
    gphit = ds["gphit"].isel(time=0, z_a=0).values
    glamu = ds["glamu"].isel(time=0, z_a=0).values
    gphiu = ds["gphiu"].isel(time=0, z_a=0).values
    glamv = ds["glamv"].isel(time=0, z_a=0).values
    gphiv = ds["gphiv"].isel(time=0, z_a=0).values
    # F-grid coordinates
    glamf = ds["glamf"].isel(time=0, z_a=0).values
    gphif = ds["gphif"].isel(time=0, z_a=0).values
    ds.close()
    # water-column depth
    bathy = np.sum(np.where(tmask, e3t, 0.0), axis=0)
    bathy = bathy.astype(float)
    bathy[np.sum(tmask, axis=0) == 0] = np.nan #Land = NaN
    #
    return (bathy, tmask, glamt, gphit, glamu, gphiu, glamv, gphiv, glamf, gphif)


# NEVER USED!
def is_deep(bathy, target):
    """
    True where the bathymetry is at least the target depth.
    Land is treated as shallow.
    """
    isdeep = (np.isfinite(bathy) & (bathy >= target))
    return isdeep


def face_crosses(a, b, target):
    """
    Return True if an edge separates a shallow point from
    a deep point.
    NaN/land is treated as shallow.
    This is important at the coastline: a deep ocean cell
    next to land must not cause the contour algorithm to
    simply disappear at the land boundary.
    """
    a_deep = (np.isfinite(a) and a >= target)
    b_deep = (np.isfinite(b) and b >= target)
    #
    return a_deep != b_deep


def build_isobath_edges(bathy, target):
    """
    Find all U and V grid faces crossed by the isobath.
    NEMO C-grid convention:
        U(j,i) is an east/west T-grid interface
        V(j,i) is a north/south T-grid interface
        U(j,i) separates T(j,i) and T(j,i+1)
        V(j,i) separates T(j,i) and T(j+1,i)
    Returns
    -------
    u_edges : bool, shape (ny,nx)
        True where U(j,i) belongs to the isobath.
    v_edges : bool, shape (ny,nx)
        True where V(j,i) belongs to the isobath.
    """
    ny, nx = bathy.shape
    u_edges = np.zeros((ny, nx), dtype=bool)
    v_edges = np.zeros((ny, nx), dtype=bool)
    # U points
    for j in range(ny):
        for i in range(nx-1):
            if face_crosses(bathy[j, i], bathy[j, i+1], target):
                u_edges[j, i] = True
    # V points
    for j in range(ny-1):
        for i in range(nx):
            if face_crosses(bathy[j, i], bathy[j+1, i], target):
                v_edges[j, i] = True
    #
    return u_edges, v_edges


def build_graph(u_edges, v_edges):
    """
    Build a graph whose vertices are the actual U/V grid edges.
    Two contour faces are connected if they meet at the same
    F-grid corner.
    This is deliberately a graph of GRID EDGES, rather than
    a graph based on longitude/latitude proximity.
    """
    ny, nx = u_edges.shape
    graph = {}
    #
    def add_node(node):
        if node not in graph:
            graph[node] = set()
    #
    def connect(a, b):
        add_node(a)
        add_node(b)
        graph[a].add(b)
        graph[b].add(a)

    # ------------------------------------------------------------
    # Each U edge:
    #     F(j,i)
    #       |
    #       |
    #     F(j-1,i)
    # ------------------------------------------------------------
    for j in range(1, ny):
        for i in range(nx):
            if u_edges[j, i]:
                f0 = ("F", j-1, i)
                f1 = ("F", j, i)
                connect(f0, f1)

    # ------------------------------------------------------------
    # Each V edge:
    #     F(j,i-1) -------- F(j,i)
    # ------------------------------------------------------------
    for j in range(ny):
        for i in range(1, nx):
            if v_edges[j, i]:
                f0 = ("F", j, i-1)
                f1 = ("F", j, i)
                connect(f0, f1)

    return graph


def connected_components(graph):
    """
    Find connected components of the contour-edge graph.
    """
    visited = set()
    components = []
    for start in graph:
        if start in visited:
            continue
        stack = [start]
        component = set()

        while stack:
            node = stack.pop()
            if node in visited:
                continue
            visited.add(node)
            component.add(node)

            for neighbour in graph[node]:
                if neighbour not in visited:
                    stack.append(neighbour)

        components.append(component)

    return components


def component_edges(component, u_edges, v_edges):
    """
    Find all U/V contour edges belonging to a connected
    component.
    """
    ny, nx = u_edges.shape
    edges = []

    # U edges
    for j in range(1, ny):
        for i in range(nx):
            if u_edges[j, i]:
                a = ("F", j-1, i)
                b = ("F", j, i)
                if a in component and b in component:
                    edges.append(("U", j, i))

    # V edges
    for j in range(ny):
        for i in range(1, nx):
            if v_edges[j, i]:
                a = ("F", j, i-1)
                b = ("F", j, i)
                if a in component and b in component:
                    edges.append(("V", j, i))

    return edges


def order_edges(edges):
    """
    Order the U/V grid edges into a continuous sequence.
    Edges are connected through common F-grid corners.
    This is only used to provide a useful ordering in the
    output NetCDF. The geometry itself is defined by the
    U/V edge indices.
    """
    if not edges:
        return []

    # ------------------------------------------------------------
    # Map each F-grid corner to contour edges touching it.
    # ------------------------------------------------------------

    endpoint_map = {}

    def add_endpoint(endpoint, edge_index):
        if endpoint not in endpoint_map:
            endpoint_map[endpoint] = []
        endpoint_map[endpoint].append(edge_index)

    for n, (grid, j, i) in enumerate(edges):

        if grid == "U":
            a = ("F", j-1, i)
            b = ("F", j, i)

        else:
            a = ("F", j, i-1)
            b = ("F", j, i)

        add_endpoint(a, n)
        add_endpoint(b, n)

    # ------------------------------------------------------------
    # Find a starting edge.
    #
    # For a closed loop every endpoint has degree 2.
    # For an open contour, start at degree 1.
    # ------------------------------------------------------------

    start = None

    for endpoint, edge_indices in endpoint_map.items():

        if len(edge_indices) == 1:
            start = edge_indices[0]
            break

    if start is None:
        start = 0

    ordered = []
    used = set()
    current = start
    while current not in used:

        used.add(current)
        ordered.append(edges[current])
        grid, j, i = edges[current]

        if grid == "U":
            endpoints = [("F", j-1, i), ("F", j, i)]
        else:
            endpoints = [("F", j, i-1), ("F", j, i)]

        # --------------------------------------------------------
        # Find the endpoint at which we should continue.
        # --------------------------------------------------------

        next_edge = None

        for endpoint in endpoints:

            for candidate in endpoint_map[endpoint]:

                if candidate != current:
                    if candidate not in used:
                        next_edge = candidate
                        break

            if next_edge is not None:
                break

        if next_edge is None:
            break

        current = next_edge

    return ordered


def write_output(filename, contours, glamu, gphiu, glamv, gphiv, target):
    """
    Write the ordered U/V isobath edges.
    """
    point_grid = []
    point_i = []
    point_j = []
    point_lon = []
    point_lat = []

    loop_start = []
    loop_length = []
    loop_closed = []

    offset = 0

    for contour in contours:

        loop_start.append(offset)
        loop_length.append(len(contour))

        if len(contour) > 1:
            closed = (contour[0] == contour[-1])
        else:
            closed = False

        loop_closed.append(closed)

        for grid, j, i in contour:

            if grid == "U":
                lon = glamu[j, i]
                lat = gphiu[j, i]
                grid_type = 0
            else:
                lon = glamv[j, i]
                lat = gphiv[j, i]
                grid_type = 1

            point_grid.append(grid_type)
            point_i.append(i)
            point_j.append(j)
            point_lon.append(lon)
            point_lat.append(lat)

            offset += 1

    ds = xr.Dataset(
        data_vars={
            "grid_type": ("point", np.asarray(point_grid, dtype=np.int8,),
                          {"long_name": ("C-grid edge type"),
                           "flag_values": [0, 1],
                           "flag_meanings": "U V",},),
            "i": ("point", np.asarray(point_i, dtype=np.int32,),),
            "j": ("point", np.asarray(point_j, dtype=np.int32,),),
            "lon": ("point", np.asarray(point_lon),),
            "lat": ("point", np.asarray(point_lat),),
            "loop_start": ("loop", np.asarray(loop_start, dtype=np.int32,),),
            "loop_length": ("loop", np.asarray(loop_length, dtype=np.int32,),),
            "closed": ("loop", np.asarray(loop_closed, dtype=bool,),),
        },

        coords={
            "point": np.arange(len(point_grid), dtype=np.int32,),
            "loop": np.arange(len(contours), dtype=np.int32,),
        },

        attrs={
            "title": "NEMO C-grid isobath",
            "isobath_depth_m": float(target),
            "description": ("Isobath represented by actual NEMO U and V grid edges."),
            "grid_type_0": "U edge",
            "grid_type_1": "V edge",
            "indexing": "zero-based",
        },
    )

    ds.to_netcdf(
        filename,
        encoding={
            "grid_type": {"dtype": "int8"},
            "i": {"dtype": "int32"},
            "j": {"dtype": "int32"},
            "loop_start": {"dtype": "int32"},
            "loop_length": {"dtype": "int32"},
        },
    )

    #ds.close()
    return ds


def compute_isobath(meshmask_file, depth, verbose=False):
    """
    Full isobath-extraction pipeline, usable as a library call.
    Reads the meshmask, finds the U/V grid faces crossed by the
    target depth and groups them into ordered contours.

    Returns a dict with the "contours" (as produced by order_edges,
    one list of (grid, j, i) edges per connected component) and the
    mesh arrays needed downstream (indexing into Fluxes.nc, plotting).
    """
    (bathy, tmask, glamt, gphit, glamu, gphiu, glamv, gphiv, glamf, gphif) = read_meshmask(meshmask_file)

    if verbose:
        print("Bathymetry range: {:.2f} - {:.2f} m".format(np.nanmin(bathy), np.nanmax(bathy)))
        print("Target depth: {:.2f} m".format(depth))

    u_edges, v_edges = build_isobath_edges(bathy, depth)

    if verbose:
        print("U isobath edges: {}".format(int(u_edges.sum())))
        print("V isobath edges: {}".format(int(v_edges.sum())))

    graph = build_graph(u_edges, v_edges)
    components = connected_components(graph)

    if verbose:
        print("Connected components: {}".format(len(components)))

    contours = []
    for component in components:
        edges = component_edges(component, u_edges, v_edges)
        ordered = order_edges(edges)
        if ordered:
            contours.append(ordered)

    if verbose:
        print("Output contours: {}".format(len(contours)))
        for n, contour in enumerate(contours):
            n_u = sum(edge[0] == "U" for edge in contour)
            n_v = sum(edge[0] == "V" for edge in contour)
            print("  {:4d}: {:6d} edges (U={:6d}, V={:6d})".format(n, len(contour), n_u, n_v))

    return {
        "contours": contours,
        "bathy": bathy,
        "glamt": glamt, "gphit": gphit,
        "glamu": glamu, "gphiu": gphiu,
        "glamv": glamv, "gphiv": gphiv,
        "glamf": glamf, "gphif": gphif,
    }


def contours_to_points(contours):
    """
    Flatten ordered U/V contours into flat point arrays.
    grid_type: 0=U edge, 1=V edge (matches write_output's convention).
    """
    point_grid = []
    point_i = []
    point_j = []

    for contour in contours:
        for grid, j, i in contour:
            point_grid.append(0 if grid == "U" else 1)
            point_i.append(i)
            point_j.append(j)

    return (np.asarray(point_i, dtype=np.int32),
            np.asarray(point_j, dtype=np.int32),
            np.asarray(point_grid, dtype=np.int8))


def main():

    parser = argparse.ArgumentParser(description=("Extract a NEMO C-grid isobath from a meshmask."))

    parser.add_argument("-m", "--meshmask")
    parser.add_argument("-d", "--depth", type=float, default=200.0)
    parser.add_argument("-o", "--outfile", default="./isobath_.nc")
    args = parser.parse_args()

    isobath = compute_isobath(args.meshmask, args.depth, verbose=True)

    # ============================================================
    # Write
    # ============================================================

    ds = write_output(args.outfile, isobath["contours"], isobath["glamu"], isobath["gphiu"],
                       isobath["glamv"], isobath["gphiv"], args.depth)

    print()
    print("Wrote {}".format(args.outfile))

    return ds


if __name__ == "__main__":
    main()

