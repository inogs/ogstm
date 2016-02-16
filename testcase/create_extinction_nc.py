import os,sys

import numpy as np

from mydtype import *

import scipy.io.netcdf as NC

import pickle

def create_extinction_nc(test):
    indata   = np.loadtxt('KB/KextF.dat', dtype=ext_data)

    jpi=test['jpi']
    jpj=test['jpj']
    
    maskfile=test['Dir'] + '/meshmask.nc'
    M=NC.netcdf_file(maskfile,"r")
    Lon     =  M.variables['glamt'].data[0,0,:,:].copy()
    Lat     =  M.variables['gphit'].data[0,0,:,:].copy()
    mask2D  =  M.variables['tmask'].data[0,0,:,:].copy()
    M.close()

    kextfact=np.ones((jpj,jpi),np.float);

    os.system("mkdir -p " + test['Dir'] + '/FORCINGS')

    for dd in indata:
        
        outfile = test['Dir'] + '/FORCINGS/KextF_' + dd['date'] + '.nc'
#        print 'creating' + outfile
        ncOUT   = NC.netcdf_file(outfile,"w")

        ncOUT.createDimension('x',jpi);
        ncOUT.createDimension('y',jpj);

        kextfact[:,:] = dd['kext'];

        ncvar = ncOUT.createVariable('nav_lon' ,'f',('y','x')); ncvar[:] = Lon     ;
        ncvar = ncOUT.createVariable('nav_lat' ,'f',('y','x')); ncvar[:] = Lat     ;
        ncvar = ncOUT.createVariable('kextfact','f',('y','x')); ncvar[:] = kextfact;

        setattr(ncOUT.variables['kextfact'],'missing_value',1e+20                    );
        setattr(ncOUT                      ,'Author'       ,'Echo group'             );
        setattr(ncOUT                      ,'Mesh'         ,'TEST_CASE'              );
        setattr(ncOUT                      ,'Origin'       ,'GOS'                    ); 
        setattr(ncOUT                      ,'Description'  ,'Extinction factor+30.0%');
        setattr(ncOUT                      ,'Units'        ,'m^(-1)'                 );
   
        ncOUT.close()
