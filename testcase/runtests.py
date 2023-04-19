import argparse
import os
from pathlib import Path
import netCDF4 as nc
import numpy as np

parser = argparse.ArgumentParser(prog=__file__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('dir1', type=Path, help="Path of first dir")
parser.add_argument('dir2', type=Path, help="Path of second dir")
parser.add_argument('--rtol', help="Tolerance used in comparison", default=1e-5, type=float)
parser.add_argument('--atol', help="Tolerance used in comparison", default=1e-8, type=float)
parser.add_argument('--log', help="Log level", default='warning')
parser.add_argument('--progress', help="Show progress bar", default=False, action='store_true')

cli_args = parser.parse_args()

def getfiles(directory):
    return {direntry.name: direntry.path 
            for direntry in filter(lambda de: de.name.endswith('.nc') and de.is_file(), os.scandir(directory))}

def intersect(A, B):
    return set.intersection(set(A), set(B))

def allclose(xs, ys):
    return np.allclose(xs, ys, rtol=cli_args.rtol, atol=cli_args.atol)

def compare_dumps(subdirectory):
    files_1, files_2 = [getfiles(directory / subdirectory) for directory in [cli_args.dir1, cli_args.dir2]]
    common_files = intersect(files_1.keys(), files_2.keys())
    tests = {}
    for filename in common_files:
        ds_1, ds_2 = [nc.Dataset(files[filename], 'r') for files in [files_1, files_2]]
        variables = filter(lambda v: v not in {'lat', 'lon', 'depth'}, intersect(ds_1.variables, ds_2.variables))
        for varname in variables:
            var_1, var_2 = [ds[varname][:] for ds in [ds_1, ds_2]]
            ispassed = allclose(var_1, var_2)
            tests[varname] = ispassed
            err = np.abs(var_1 - var_2).max()
            print(f"{varname} ({err=:.2e}):\t{'PASS' if ispassed else 'FAIL'}")

    print(32 * "=")
    if all(tests.values()):
        print(f"All tests passed for {subdirectory}")
    else:
        print(f"Some tests failed for {subdirectory}")
    print(32 * "=")

if __name__ == '__main__':
    for subdir in ['AVE_FREQ_1', 'AVE_FREQ_2']:
        compare_dumps(subdir)


