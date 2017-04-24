import os,sys

import numpy as np

from mydtype import *

import scipy.io.netcdf as NC

import pickle

import datetime

from scipy.interpolate import interp1d

def create_TSKQWHF(test,flnm,D3T,D3SIGMA,D3S,D3K,D2W,D2H,D2F,PAR,Ed380,Ed412,Ed490,CHL_F):

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

    float_data=np.loadtxt(flnm, dtype=Bio_Float_type,skiprows=1,ndmin=1)

#   getting Julian date for present simulation

    gdeptTOT= float_data['Depth']
    s_0     = float_data['SAL'];# Input data to be interpolated on final grid
    t_0     = float_data['Temp'];# Input data to be interpolated on final grid
    si_0    = float_data['sigma'];# Input data to be interpolated on final grid
    k_0     = float_data['EKD'];# Input data to be interpolated on final grid
    print('PAR is divided by a factor 2.')
    par_0   = float_data['PAR']*0.5;
    Ed380_0 = float_data['Ed380'];
    Ed412_0 = float_data['Ed412'];
    Ed490_0 = float_data['Ed490'];
    CHL_F_0 = float_data['CHL_F'];

#   rho   = 1.3 # kg/m3
#   Cdrag = 1.5 * 0.001
#   K     = np.sqrt(1./(rho*Cdrag))
#   w     = (taux**2 + tauy**2)**0.25 * K

    w  =0.5
    h  =0.;
    f  =0.;
#   f  =DATA_T.variables['sowaflcd'].data[j_day,0,0].copy();

# Interpolating from original vertical grid to high resolution one

    t_int     = interp1d(gdeptTOT,t_0,kind='nearest',fill_value='extrapolate');     t = t_int(gdept)
    s_int     = interp1d(gdeptTOT,s_0,kind='nearest',fill_value='extrapolate');     s = s_int(gdept) 
    si_int    = interp1d(gdeptTOT,si_0,kind='nearest',fill_value='extrapolate');    si = si_int(gdept) 
    k_int     = interp1d(gdeptTOT,k_0,kind='nearest',fill_value='extrapolate');     k = k_int(gdepw)
    par_int   = interp1d(gdeptTOT,par_0,kind='nearest',fill_value='extrapolate');   par_1d = par_int(gdept)
    Ed380_int = interp1d(gdeptTOT,Ed380_0,kind='nearest',fill_value='extrapolate'); Ed380_1d = Ed380_int(gdept) 
    Ed412_int = interp1d(gdeptTOT,Ed412_0,kind='nearest',fill_value='extrapolate'); Ed412_1d = Ed412_int(gdept)
    Ed490_int = interp1d(gdeptTOT,Ed490_0,kind='nearest',fill_value='extrapolate'); Ed490_1d = Ed490_int(gdept)
    CHL_F_int = interp1d(gdeptTOT,CHL_F_0,kind='nearest',fill_value='extrapolate'); CHL_F_1d = CHL_F_int(gdept)

##############################

    for jk in np.arange(jpk):
       for jj in np.arange(jpj):
           for ji in np.arange(jpi):
               D3S[0,jk,jj,ji] = s[jk]
               D3T[0,jk,jj,ji] = t[jk]
               D3SIGMA[0,jk,jj,ji] = si[jk]
               if jk > 4 :
                   D3K[0,jk,jj,ji] = np.amin([k[jk],D3K[0,jk-1,jj,ji]])
                   D3K[0,jk,jj,ji] = np.amax([D3K[0,jk,jj,ji],0.000001])
               else:
                   D3K[0,jk,jj,ji] = np.amax([k[jk],0.000001])
               PAR[0,jk,jj,ji]   = par_1d[jk]
               Ed380[0,jk,jj,ji] = Ed380_1d[jk]
               Ed412[0,jk,jj,ji] = Ed412_1d[jk]
               Ed490[0,jk,jj,ji] = Ed490_1d[jk]
               CHL_F[0,jk,jj,ji] = CHL_F_1d[jk]

    for jj in np.arange(jpj):
        for ji in np.arange(jpi):
            D2W[0,jj,ji] = w
            D2H[0,jj,ji] = 0.
            D2F[0,jj,ji] = f
