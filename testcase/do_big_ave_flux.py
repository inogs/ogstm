import os,sys,glob

import numpy as np

from mydtype import *

import scipy.io.netcdf as NC

import pickle

import datetime

def do_big_ave_flux(test):

    jpi=test['jpi'];
    jpj=test['jpj'];
    jpk=test['jpk'];
    time = 1
    maskfile=test['Dir'] + '/meshmask.nc'

    M=NC.netcdf_file(maskfile,"r")

    Lon     =  M.variables['glamt'].data[0,0,:,:].copy()
    Lat     =  M.variables['gphit'].data[0,0,:,:].copy()
    gdept   =  M.variables['gdept'].data[0,:,0,0].copy()
    gdepw   =  M.variables['gdepw'].data[0,:,0,0].copy()

    M.close()

# Retrieving state variables names from model namelist
    CODEPATH = test['Code'] + "/ogstm/"
    CODEPATH = CODEPATH.replace("~",os.getenv("HOME"))
    filename = CODEPATH +  "ready_for_model_namelists/namelist.passivetrc"
    NAMELIST = file2stringlist(filename)
    VARS=[]
    for line in NAMELIST:
         if line.find("ctrcnm") != -1:
            quote_1=line.find("\"")
            quote_2=line.find("\"",quote_1+1)
            varname=line[quote_1+1:quote_2]
            VARS.append(varname)

    nvars=len(VARS)
    check_bool = 0
    for v,var in enumerate(VARS):
         DIR_DATA          = test['Dir'] + '/'
         print DIR_DATA 
#/pico/scratch/userexternal/plazzari/TILMAN/ogstm/testcase
#        DIR_DATA          = '/pico/scratch/userexternal/plazzari/TILMAN/ogstm/testcase/TEST01/wrkdir/MODEL/' + 'AVE_FREQ_1/'
         filenames         = DIR_DATA + 'flux.2002*.nc'
         SingleVar_filelist= glob.glob(filenames)
         SingleVar_filelist.sort()
         if check_bool == 0:
             ntimes=len(SingleVar_filelist)
             matrix=np.zeros((ntimes,jpk-1,7))
             #       WRITE NetCDF restart file
             outfile = 'POSTPROC/' + test['Area'] + '_flux.nc'
             ncOUT   = NC.netcdf_file(outfile,"w");
        
             ncOUT.createDimension('depth',jpk-1);
             ncOUT.createDimension('time',ntimes)
             ncOUT.createDimension('Type',7)
             check_bool = 1
         for nf,fname in enumerate(SingleVar_filelist):
             DATA=NC.netcdf_file(fname,"r")
             matrix[nf,:,:]=DATA.variables[var].data[:,:].copy();# Input data to be interpolated on final grid
             DATA.close()
         ncvar = ncOUT.createVariable(var       ,'f',('time','depth','Type')); ncvar[:] = matrix; 
         setattr(ncOUT.variables[var]   ,'missing_value',1e+20                              );     
    ncOUT.close()
