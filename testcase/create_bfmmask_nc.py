#! /usr/bin/python

# LOAD PACKAGES
import os,sys

import numpy as np

from mydtype import *

import scipy.io.netcdf as NC

import pickle

# Script to create bfm mask file


def create_bfm_nc(test):
    x=test['jpi']
    y=test['jpj']
    z=test['jpk']

    mask = np.ones((z,y,x),bool);

    ##############################################################
    # write meshmask netcdf file !
    ##############################################################
    os.system("mkdir -p " + test['Dir'].decode())

    outfile = test['Dir'].decode() + '/bfmmask.nc';

    ncOUT=NC.netcdf_file(outfile,"w");

    ncOUT.createDimension('x',x);
    ncOUT.createDimension('y',y);
    ncOUT.createDimension('z',z);

    ncvar    = ncOUT.createVariable('bfmmask','b',('z','y','x'))  ; ncvar[:] = mask;
    ncOUT.close()

