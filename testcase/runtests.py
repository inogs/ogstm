#! /usr/bin/env python

import argparse
import os
from pathlib import Path
from tqdm import tqdm
from colorama import Fore as fore
from colorama import Style as stl
from colorama import init as colorama_init
import netCDF4 as nc
import numpy as np

parser = argparse.ArgumentParser(prog=__file__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('dir1', type=Path, help="Path of first dir")
parser.add_argument('dir2', type=Path, help="Path of second dir")
parser.add_argument('subdir', type=Path, help="Path of common subdir")
parser.add_argument('--rtol', help="Tolerance used in comparison", default=1e-5, type=float)
parser.add_argument('--atol', help="Tolerance used in comparison", default=1e-8, type=float)
parser.add_argument('--noprogress', help="Disable progress bar", default=False, action='store_true')

cli_args = parser.parse_args()

def getfiles(directory):
    return {direntry.name: direntry.path 
            for direntry in filter(lambda de: de.name.endswith('.nc') and de.is_file(), os.scandir(directory))}

def intersect(A, B):
    return set.intersection(set(A), set(B))

def allclose(xs, ys):
    return np.ma.allclose(xs, ys, rtol=cli_args.rtol, atol=cli_args.atol)

def colorize(string, passed=True):
    return stl.BRIGHT + (fore.GREEN if passed else fore.RED) + string + stl.RESET_ALL

def compare_dumps(subdirectory):
    files_1, files_2 = [getfiles(directory / subdirectory) for directory in [cli_args.dir1, cli_args.dir2]]
    common_files = intersect(files_1.keys(), files_2.keys())
    tests = {}
    for filename in tqdm(common_files, disable=cli_args.noprogress):
        ds_1, ds_2 = [nc.Dataset(files[filename], 'r') for files in [files_1, files_2]]
        variables = filter(lambda v: v not in {'lat', 'lon', 'depth'}, intersect(ds_1.variables, ds_2.variables))
        for varname in variables:
            var_1, var_2 = [ds[varname][:] for ds in [ds_1, ds_2]]
            ispassed = allclose(var_1, var_2)
            diffs = np.ma.abs(var_1 - var_2)
            tests[varname] = ispassed, diffs.max()
    
    return tests

def print_results(tests):
    print('\n' + 'VARIABLE' + 8 * ' ' + 'MAX DIFF' + 9 * ' ' + 'RESULT' + 26 * ' ')
    print(48 * "=")
    for varname, (ispassed, maxdiff) in tests.items():
        print(f"{varname:16}{maxdiff:.2e}" + 9 * ' ', end='')
        print(colorize('PASS') if ispassed else colorize('FAIL', passed=False))
    if all(ispassed for ispassed, _ in tests.values()):
        print(colorize(f"\nAll tests passed"))
    else:
        print(colorize(f"\nSome tests failed", passed=False))

colorama_init()

if __name__ == '__main__':
    tqdm.write("Running tests\n")
    tests = compare_dumps(cli_args.subdir)
    print_results(tests)
   

