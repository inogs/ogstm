#!/usr/bin/env python

import argparse

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

import define_isobath

'''
self explaining

Can also be imported as a library: plot_isobath_map(contours, bathy,
lon_t, lat_t, lon_f, lat_f, ...) draws the bathymetry + isobath overlay
on a given/new axis (see create_fluxes.py).

Usage: 
source /g100/home/userexternal/ggalli00/PY_ENVS/bitsea_gg/sequence3.sh
isobath_file='/g100_work/OGS_test2528/ggalli00/EIS_2026_11/SCRIPTS/isobath_200m.nc'
maskfile='/g100_work/OGS_test2528/ggalli00/EIS_2026_11/qDEG_setup/MASK/meshmask_025_z125.nc'
python plot_isobath.py $isobath_file $maskfile 
'''


def read_isobath(filename):

    ds = xr.open_dataset(filename)

    grid_type = ds["grid_type"].values
    ii = ds["i"].values
    jj = ds["j"].values

    loop_start = ds["loop_start"].values
    loop_length = ds["loop_length"].values

    contours = []

    for start, length in zip(
        loop_start,
        loop_length
    ):

        start = int(start)
        stop = start + int(length)

        contour = []

        for grid, j, i in zip(
            grid_type[start:stop],
            jj[start:stop],
            ii[start:stop]
        ):

            name = (
                "U"
                if int(grid) == 0
                else "V"
            )

            contour.append(
                (
                    name,
                    int(j),
                    int(i)
                )
            )

        contours.append(contour)

    ds.close()

    return contours


def plot_u_edge(
    ax,
    j,
    i,
    lon_f,
    lat_f,
    **kwargs
):
    """
    U[j,i] is the east face of T[j,i].
    It runs between F[j-1,i] and F[j,i].
    """

    if j <= 0:
        return

    if j >= lon_f.shape[0]:
        return

    if i < 0 or i >= lon_f.shape[1]:
        return

    ax.plot(
        [lon_f[j-1, i], lon_f[j, i]],
        [lat_f[j-1, i], lat_f[j, i]],
        **kwargs
    )


def plot_v_edge(
    ax,
    j,
    i,
    lon_f,
    lat_f,
    **kwargs
):
    """
    V[j,i] is the north face of T[j,i].
    It runs between F[j,i-1] and F[j,i].
    """

    if j < 0 or j >= lon_f.shape[0]:
        return

    if i < 0 or i + 1 >= lon_f.shape[1]:
        return

    ax.plot(
        [lon_f[j, i-1], lon_f[j, i]],
        [lat_f[j, i-1], lat_f[j, i]],
        **kwargs
    )


def plot_isobath(
    ax,
    contours,
    lon_f,
    lat_f
):

    for contour in contours:

        for grid, j, i in contour:

            if grid == "U":

                plot_u_edge(
                    ax,
                    j,
                    i,
                    lon_f,
                    lat_f,
                    color="red",
                    linewidth=2.0,
                    solid_capstyle="butt",
                    zorder=10
                )

            else:

                plot_v_edge(
                    ax,
                    j,
                    i,
                    lon_f,
                    lat_f,
                    color="red",
                    linewidth=2.0,
                    solid_capstyle="butt",
                    zorder=10
                )


def plot_isobath_map(
    contours,
    bathy,
    lon_t,
    lat_t,
    lon_f,
    lat_f,
    title="NEMO isobath",
    output=None,
    ax=None
):
    """
    Draw the bathymetry (pcolormesh) with the isobath contour(s)
    overlaid. Library entry point, usable without going through the
    isobath/meshmask netcdf files (see create_fluxes.py).

    If ax is None a new figure/axis is created, otherwise the plot is
    drawn on the given axis. If output is given, the figure is saved
    there. Returns (fig, ax).
    """

    if ax is None:
        fig, ax = plt.subplots(figsize=(11, 7))
    else:
        fig = ax.figure

    pcm = ax.pcolormesh(
        lon_t,
        lat_t,
        np.ma.masked_invalid(bathy),
        shading="auto",
        cmap="Blues_r"
    )

    cbar = fig.colorbar(
        pcm,
        ax=ax
    )

    cbar.set_label(
        "Depth (m)"
    )

    plot_isobath(
        ax,
        contours,
        lon_f,
        lat_f
    )

    ax.set_xlabel(
        "Longitude"
    )

    ax.set_ylabel(
        "Latitude"
    )

    ax.set_title(
        title
    )

    ax.set_aspect(
        "equal",
        adjustable="box"
    )

    fig.tight_layout()

    if output:
        fig.savefig(
            output,
            dpi=200,
            bbox_inches="tight"
        )

    return fig, ax


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("isobath")

    parser.add_argument("meshmask")

    parser.add_argument("-o", "--output", default=None)

    args = parser.parse_args()

    (bathy, tmask, glamt, gphit, glamu, gphiu, glamv, gphiv, glamf, gphif) = define_isobath.read_meshmask(args.meshmask)

    contours = read_isobath(args.isobath)

    fig, ax = plot_isobath_map(
        contours,
        bathy,
        glamt,
        gphit,
        glamf,
        gphif,
        title="NEMO 200 m isobath",
        output=args.output
    )

    if not args.output:
        plt.show()


if __name__ == "__main__":
    main()

