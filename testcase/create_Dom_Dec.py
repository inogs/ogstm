import os,sys

import numpy as np

from mydtype import *

import scipy.io.netcdf as NC

import pickle
import importlib.util

spec = importlib.util.spec_from_file_location('domdec', '../preproc/domdec/domdec.py')
domdec = importlib.util.module_from_spec(spec)
spec.loader.exec_module(domdec)

def create_Dom_Dec(test):

    jpi = test['jpi']
    jpj = test['jpj']

    nPx = test['nprocx']
    nPy = test['nprocy']

    nx = np.zeros(nPx,dtype='int')
    ny = np.zeros(nPy,dtype='int')

    for ji in range(nPx):
        nx[ji]  = jpi/nPx + 2; # 2 ghost cell

    nx[0] +=         -1 # first and last only one ghost cell
    nx[-1]+= jpi%nPx -1

    for jj in range(nPy):
        ny[jj]  = jpj/nPy + 2; # 2 ghost cell

    ny[0] +=         -1
    ny[-1]+= jpj%nPy -1

    os.system("mkdir -p " + test['Dir'].decode())

    filename01 = test['Dir'].decode() + '/Dom_Dec_jpi.ascii'
    filename02 = test['Dir'].decode() + '/Dom_Dec_jpj.ascii'

    f01 = open(filename01,'w')
    f02 = open(filename02,'w')

    for ji in range(nPx):
        for jj in range(nPy):

            s1 = str(nx[ji])  + ' '
            s2 = str(ny[jj])  + ' '
            f01.write(s1)
            f02.write(s2)

        f01.write("\n")
        f02.write("\n")

    f01.close()
    f02.close()
    tmask = np.ones((jpj,jpi),dtype=bool)
    domdec.dump_outfile(tmask,nPx*nPy, nPx, nPy, filename=test['Dir'].decode() + '/domdec.txt')

