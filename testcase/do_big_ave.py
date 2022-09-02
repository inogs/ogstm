import os,sys,glob

import numpy as np

from mydtype import *

import scipy.io.netcdf as NC
import netCDF4 as NC4

import pickle

import datetime

def do_big_ave(test,POSTPROC_DIR):

    jpi=test['jpi'];
    jpj=test['jpj'];
    jpk=test['jpk'];
    time = 1
    print(test['Dir'])
    maskfile=test['Dir'].decode() + "/meshmask.nc"

    M=NC.netcdf_file(maskfile,"r")

    Lon     =  M.variables['glamt'].data[0,0,:,:].copy()
    Lat     =  M.variables['gphit'].data[0,0,:,:].copy()
    gdept   =  M.variables['gdept'].data[0,:,0,0].copy()
    gdepw   =  M.variables['gdepw'].data[0,:,0,0].copy()

    M.close()

# Retrieving state variables names from model namelist
    CODEPATH = test['Code'].decode() + "/ogstm/"
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
         DIR_DATA          = test['Dir'].decode() + '/AVE_FREQ_2/'
         print(DIR_DATA)
#/pico/scratch/userexternal/plazzari/TILMAN/ogstm/testcase
#        DIR_DATA          = '/pico/scratch/userexternal/plazzari/TILMAN/ogstm/testcase/TEST01/wrkdir/MODEL/' + 'AVE_FREQ_1/'
         filenames         = DIR_DATA + 'ave.*.' + var +'.nc'
         SingleVar_filelist= glob.glob(filenames)
         SingleVar_filelist.sort()
         if check_bool == 0:
             ntimes=len(SingleVar_filelist)
             matrix=np.zeros((ntimes,jpk))
             #       WRITE NetCDF restart file
             # outfile = POSTPROC_DIR + '/' + test['BIO-FLOAT'].decode() + '.nc'
             outfile = POSTPROC_DIR + '/' + 'TEST01' + '.nc'
             ncOUT   = NC4.Dataset(outfile,"w");
        
             ncOUT.createDimension('depth',jpk);
             ncOUT.createDimension('time',ntimes);
             check_bool = 1
         for nf,fname in enumerate(SingleVar_filelist):
             DATA=NC4.Dataset(fname,"r")
#            DATA=NC.netcdf_file(fname,"r")
             matrix[nf,:]=DATA.variables[var][0,:,1,1].copy();# Input data to be interpolated on final grid
             DATA.close()
         ncvar = ncOUT.createVariable(var       ,'f',('time','depth'),zlib=True); ncvar[:] = matrix; 
         setattr(ncOUT.variables[var]   ,'missing_value',1e+20                              );     
    ncOUT.close()


#   getting Julian date for present simulation
#   mm    = date[4:6]
#   dd    = date[6:8]
#   s0    = mm+'.'+dd
#   fmt   = '%m.%d'
#   dt    = datetime.datetime.strptime(s0,fmt)
#   tt    = dt.timetuple()
#   j_day = tt.tm_yday -1 # from 1 to 365 scaled to 0 364

