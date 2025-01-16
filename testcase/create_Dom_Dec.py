import os,sys

import numpy as np

from mydtype import *
from importlib.machinery import SourceFileLoader
domdec = SourceFileLoader("domdec", "../preproc/domdec/domdec.py").load_module()


def create_Dom_Dec(test):

    jpi = test['jpi']
    jpj = test['jpj']

    nPx = test['nprocx']
    nPy = test['nprocy']

    os.system("mkdir -p " + test['Dir'])


    tmask = np.ones((jpj,jpi),dtype=bool)
    domdec.dump_outfile(tmask,nPx*nPy, nPx, nPy, filename=test['Dir'] + '/domdec.txt')


