import os,sys

import numpy as np

from mydtype import *

import scipy.io.netcdf as NC

import pickle


def h_vort1(jpi,jpj,jpk,D3U,D3V,Av,Au):

    if (jpj <= jpi):
        if jpj%2 == 0:
            ncircle=(jpj-2)/2
        else:
            ncircle=(jpj-3)/2
    else:
        if jpi%2 == 0:
            ncircle=(jpi-2)/2
        else:
            ncircle=(jpi-3)/2
# tangent flux
    flux=0.5 *10**6; # Flux in sverdrup m^3/s
    for n in range(ncircle):
        
        for jk in range(jpk):
        # Velocities along jpj side Ly
            Ly = np.arange(1+n,jpj-2-n , dtype=np.int)
            for jj in Ly:
                D3V[0,jk,jj,1+n    ] =  flux/Av[jk,jj,1+n    ]     # one side  
                D3V[0,jk,jj,jpi-2-n] = -flux/Av[jk,jj,jpi-2-n]
            # Velocities along jpi side Lx
            Lx = np.arange(1+n,jpi-2-n , dtype=np.int)
            for ji in Lx:
                D3U[0,jk,1+n,ji]     = -flux/Au[jk,1+n,ji]     # one side 
                D3U[0,jk,jpj-2-n,ji] =  flux/Au[jk,jpj-2-n,ji]


def create_UVW(test,D3U,D3V,D3W,Av,Au):

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

    h_vort1(jpi,jpj,jpk,D3U,D3V,Av,Au)


