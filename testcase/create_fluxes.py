import numpy as np

from mydtype import *

import scipy.io.netcdf as NC

def create_fluxes(test):
    '''
    Creates Fluxes.nc from a simple transect defined inside the function
    The transect is defined as 
     - j=5
     - k<5
    '''

    jpi=test['jpi'];
    jpj=test['jpj'];
    jpk=test['jpk'];
    modeldir=test['Dir'].decode()
    
    maskfile=modeldir + '/meshmask.nc'
    bounmask=modeldir + '/bounmask.nc'
    outfile = modeldir + "/Fluxes.nc"

    M=NC.netcdf_file(maskfile,"r")
    tmask     =  M.variables['tmask'].data[0,:,:,:].copy()
    M.close()
    
    M=NC.netcdf_file(bounmask,"r")
    index     =  M.variables['index'].data[0,:,:,:].copy()
    M.close()
    
    j_transect=5
    LIST=[]
    for i in range(1,jpi-1):
        for k in range(0,5):
            new_el = index[k,j_transect,i]
            if new_el in LIST: print "already", i,k, new_el
            LIST.append(new_el) 
    
    n=len(LIST)

    ncOUT=NC.netcdf_file(outfile,"w")
    ncOUT.createDimension('n',n)
    ncvar=ncOUT.createVariable('index','i',('n',))
    ncvar[:]=LIST
    ncOUT.close()

