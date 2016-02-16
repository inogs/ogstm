import os,sys

import numpy as np

from mydtype import *

import scipy.io.netcdf as NC

import pickle


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

    forcfile= "KB/FORCINGS/" +  "ave." + date + ".stat_profiles.nc"

    DATA=NC.netcdf_file(forcfile,"r")

#   sub___list = "alb, sww, swe, nwm, tyr, adn, ads, aeg, ion, lev, med" ;
#   stat__list = "Mean, Std, p25, p50, p75" ;
#   coast_list = "coast, open_sea, everywhere" ;

#   float votemper(sub, coast, z, stat) ;

    s=DATA.variables['vosaline'].data[3,1,:,0].copy();
    t=DATA.variables['votemper'].data[3,1,:,0].copy();
    k=DATA.variables['votkeavt'].data[3,1,:,0].copy();
    q=DATA.variables['soshfldo'].data[3,1,0,0].copy();
    w=DATA.variables['sowindsp'].data[3,1,0,0].copy();
    h=0.;#DATA.variables['sosheigh'].data[3,1,0,0].copy();
    f=DATA.variables['sowaflcd'].data[3,1,0,0].copy();

# filling zero values
    s_f =s.copy()
    t_f =t.copy()
    k_f =k.copy()

    for jk in np.arange(jpk):
        if s[jk] == 0:
            idx =jk 
            break
    s_f[idx:jpk]=s[idx-1] 

    for jk in np.arange(jpk):
        if t[jk] == 0:
            idx =jk 
            break
    t_f[idx:jpk]=t[idx-1] 

    for jk in np.arange(jpk):
        if k[jk] == 0:
            idx =jk 
            break
    k_f[idx:jpk]=k[idx-1] 
      
##############################

    for jk in np.arange(jpk):
       for jj in np.arange(jpj):
           for ji in np.arange(jpi):
               D3S[0,jk,jj,ji] = s_f[jk]
               D3T[0,jk,jj,ji] = t_f[jk]
               D3K[0,jk,jj,ji] = k_f[jk]

    for jj in np.arange(jpj):
        for ji in np.arange(jpi):
            D2Q[0,jj,ji] = q
            D2W[0,jj,ji] = w
            D2H[0,jj,ji] = 0.
            D2F[0,jj,ji] = f

    DATA.close()
