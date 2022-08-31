#! /usr/bin/python

# LOAD PACKAGES
import os,sys

import numpy as np

from mydtype import *

import scipy.io.netcdf as NC

# Script to create bfmmask file

def create_bfmmask_nc(test):
    jpi=test['jpi']
    jpj=test['jpj']
    jpk=test['jpk']

    dx=test['dx']
    dy=test['dy']

#    double bfmmask( z, y, x) ;
    bfmmask = np.ones((jpk,jpj,jpi),np.bool);

    ##############################################################
    # write bfmmask netcdf file !
    ##############################################################
    os.system("mkdir -p " + test['Dir'].decode())

    outfile = test['Dir'].decode() + '/bfmmask.nc';

    ncOUT=NC.netcdf_file(outfile,"w");

    ncOUT.createDimension('x',jpi);
    ncOUT.createDimension('y',jpj);
    ncOUT.createDimension('z',jpk);

    ncvar    = ncOUT.createVariable('bfmmask'   ,'b',('z', 'y', 'x') )  ; ncvar[:] = bfmmask   ;
    ncOUT.close()

