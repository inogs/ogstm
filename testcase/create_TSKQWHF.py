import os,sys

import numpy as np

from mydtype import *

import scipy.io.netcdf as NC

import pickle

import datetime

from scipy.interpolate import interp1d

def create_TSKQWHF(test,date,D3T,D3S,D3K,D2Q,D2W,D2H,D2F):

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

    Area=test['Area']

    forcfileT= "COPERNICUS/FORCINGS/" + Area + "/CLIM/" +  Area + "_T.nc"
    forcfileW= "COPERNICUS/FORCINGS/" + Area + "/CLIM/" +  Area + "_W.nc"

    DATA_T=NC.netcdf_file(forcfileT,"r")
    DATA_W=NC.netcdf_file(forcfileW,"r")

#   getting Julian date for present simulation
    mm    = date[4:6]
    dd    = date[6:8]
    s0    = mm+'.'+dd
    fmt   = '%m.%d'
    dt    = datetime.datetime.strptime(s0,fmt)
    tt    = dt.timetuple()
    j_day = tt.tm_yday -1 # from 1 to 365 scaled to 0 364

    s_0=DATA_T.variables['vosaline'].data[j_day,:,0,0].copy();# Input data to be interpolated on final grid
    t_0=DATA_T.variables['votemper'].data[j_day,:,0,0].copy();# Input data to be interpolated on final grid
    k_0=DATA_W.variables['votkeavt'].data[j_day,:,0,0].copy();# Input data to be interpolated on final grid
    q  =DATA_T.variables['soshfldo'].data[j_day,0,0].copy();
    w  =1.
#   w  =DATA_T.variables['sowindsp'].data[j_day,0,0].copy();
    h  =0.;
    f  =0.;
#   f  =DATA_T.variables['sowaflcd'].data[j_day,0,0].copy();

# Interpolating from original vertical grid to high resolution one
    filein             = 'COPERNICUS' + '/gdept' + 'COPERNICUS' + '.dat'
    gdeptTOT           = np.loadtxt(filein, dtype=np.double);

    t_int = interp1d(gdeptTOT,t_0,kind='nearest',fill_value='extrapolate'); t = t_int(gdept)
    s_int = interp1d(gdeptTOT,s_0,kind='nearest',fill_value='extrapolate'); s = s_int(gdept) 
    
   
    filein             = 'COPERNICUS' + '/gdepw' + 'COPERNICUS' + '.dat'
    gdepwTOT           = np.loadtxt(filein, dtype=np.double);
    k_int = interp1d(gdepwTOT,k_0,kind='nearest',fill_value='extrapolate'); k = k_int(gdepw)

##############################

    for jk in np.arange(jpk):
       for jj in np.arange(jpj):
           for ji in np.arange(jpi):
               D3S[0,jk,jj,ji] = s[jk]
               D3T[0,jk,jj,ji] = t[jk]
               if jk > 4 :
                   D3K[0,jk,jj,ji] = np.amin([k[jk],D3K[0,jk-1,jj,ji]])
                   D3K[0,jk,jj,ji] = np.amax([D3K[0,jk,jj,ji],0.00000001])
               else:
                   D3K[0,jk,jj,ji] = np.amax([k[jk],0.00000001])

    for jj in np.arange(jpj):
        for ji in np.arange(jpi):
            D2Q[0,jj,ji] = q
            D2W[0,jj,ji] = w
            D2H[0,jj,ji] = 0.
            D2F[0,jj,ji] = f

    DATA_T.close()
    DATA_W.close()
