import os,sys

import numpy as np

from mydtype import *

import scipy.io.netcdf as NC

import pickle

import create_UVW as UVW
import create_TSKQWHF  as TSKQWHF

def create_forcings_nc(test):

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
    e1v     =  M.variables['e1v'  ].data[0,0,:,:].copy()
    e2u     =  M.variables['e2u'  ].data[0,0,:,:].copy()
    e3t     =  M.variables['e3t'  ].data[0,:,0,0].copy()
    
    M.close()
    
    e1v_aux = np.tile(e1v,[jpk,1,1])
    e2u_aux = np.tile(e2u,[jpk,1,1])

    Av=np.ones((jpk,jpj,jpi),np.float)
    Au=np.ones((jpk,jpj,jpi),np.float)
    
    for i in range(jpk):
        Av[i,:,:] = e1v_aux[i,:,:] * e3t[i]
        Au[i,:,:] = e2u_aux[i,:,:] * e3t[i]

   
    D3T=np.zeros((1,jpk,jpj,jpi),np.float)
    D3S=np.zeros((1,jpk,jpj,jpi),np.float)
    D3K=np.zeros((1,jpk,jpj,jpi),np.float)

    D2Q=np.zeros((1,jpj,jpi),np.float)
    D2W=np.zeros((1,jpj,jpi),np.float)
    D2H=np.zeros((1,jpj,jpi),np.float)
    D2F=np.zeros((1,jpj,jpi),np.float)

    D3U=np.zeros((1,jpk,jpj,jpi),np.float)
    D3V=np.zeros((1,jpk,jpj,jpi),np.float)
    D3W=np.zeros((1,jpk,jpj,jpi),np.float)
    
    UVW.create_UVW(test,D3U,D3V,D3W,Av,Au)
    
    D3=np.ones((1,jpk,jpj,jpi),np.float)    
    D2=np.ones((1,jpj,jpi),np.float)


    FORCING_DATE=[]

    filein=file('KB/forcing_date')
    for var in filein:
        FORCING_DATE.append(var[:-1])
    
    filein.close()    

    os.system("mkdir -p " + test['Dir'] + '/FORCINGS/')

    for date in FORCING_DATE:
        # Create T file
        TSKQWHF.create_TSKQWHF(test,date,D3T,D3S,D3K,D2Q,D2W,D2H,D2F)
        outfile = test['Dir'] + '/FORCINGS/T' + date + '.nc'
        ncOUT   = NC.netcdf_file(outfile,'w')

        ncOUT.createDimension('x'           ,jpi);
        ncOUT.createDimension('y'           ,jpj);
        ncOUT.createDimension('deptht'      ,jpk);
        ncOUT.createDimension('time_counter',time);
        
        ncvar = ncOUT.createVariable('nav_lon'      ,'f',('y','x')                        ); ncvar[:] = Lon;
        ncvar = ncOUT.createVariable('nav_lat'      ,'f',('y','x')                        ); ncvar[:] = Lat;
        ncvar = ncOUT.createVariable('deptht'       ,'f',('deptht',)                      ); ncvar[:] = gdept;
        ncvar = ncOUT.createVariable('time_counter' ,'d',('time_counter',)                ); ncvar    = 1.;
        ncvar = ncOUT.createVariable('vosaline'     ,'f',('time_counter','deptht','y','x')); ncvar[:] = D3S;   
        ncvar = ncOUT.createVariable('votemper'     ,'f',('time_counter','deptht','y','x')); ncvar[:] = D3T;
        ncvar = ncOUT.createVariable('soshfldo'     ,'f',('time_counter','y','x')         ); ncvar[:] = D2Q;
        ncvar = ncOUT.createVariable('sowindsp'     ,'f',('time_counter','y','x')         ); ncvar[:] = D2W;
        ncvar = ncOUT.createVariable('sossheig'     ,'f',('time_counter','y','x')         ); ncvar[:] = 0.;
        ncvar = ncOUT.createVariable('sowaflcd'     ,'f',('time_counter','y','x')         ); ncvar[:] = D2F;
        ncOUT.close()

        # Create U file

        outfile = test['Dir'] + '/FORCINGS/U' + date + '.nc'
        ncOUT   = NC.netcdf_file(outfile,'w')

        ncOUT.createDimension('x'           ,jpi);
        ncOUT.createDimension('y'           ,jpj);
        ncOUT.createDimension('deptht'      ,jpk);
        ncOUT.createDimension('time_counter',time);
        
        ncvar = ncOUT.createVariable('nav_lon'      ,'f',('y','x')                        ); ncvar[:] = Lon;
        ncvar = ncOUT.createVariable('nav_lat'      ,'f',('y','x')                        ); ncvar[:] = Lat;
        ncvar = ncOUT.createVariable('deptht'       ,'f',('deptht',)                      ); ncvar[:] = gdept;
        ncvar = ncOUT.createVariable('time_counter' ,'d',('time_counter',)                ); ncvar    = 1.;
        ncvar = ncOUT.createVariable('vozocrtx'     ,'f',('time_counter','deptht','y','x')); ncvar[:] = D3U; 
        ncOUT.close()

        # Create V file

        outfile = test['Dir'] + '/FORCINGS/V' + date + '.nc'
        ncOUT   = NC.netcdf_file(outfile,'w')

        ncOUT.createDimension('x'           ,jpi);
        ncOUT.createDimension('y'           ,jpj);
        ncOUT.createDimension('deptht'      ,jpk);
        ncOUT.createDimension('time_counter',time);
        
        ncvar = ncOUT.createVariable('nav_lon'      ,'f',('y','x')                        ); ncvar[:] = Lon;
        ncvar = ncOUT.createVariable('nav_lat'      ,'f',('y','x')                        ); ncvar[:] = Lat;
        ncvar = ncOUT.createVariable('deptht'       ,'f',('deptht',)                      ); ncvar[:] = gdept;
        ncvar = ncOUT.createVariable('time_counter'  ,'d',('time_counter',)               ); ncvar    = 1.;
        ncvar = ncOUT.createVariable('vomecrty'     ,'f',('time_counter','deptht','y','x')); ncvar[:] = D3V; 
        ncOUT.close()

        # Create W file

        outfile = test['Dir'] + '/FORCINGS/W' + date + '.nc'
        ncOUT   = NC.netcdf_file(outfile,'w')

        ncOUT.createDimension('x'           ,jpi);
        ncOUT.createDimension('y'           ,jpj);
        ncOUT.createDimension('depthw'      ,jpk);
        ncOUT.createDimension('time_counter',time);
        
        ncvar = ncOUT.createVariable('nav_lon'      ,'f',('y','x')                        ); ncvar[:] = Lon;
        ncvar = ncOUT.createVariable('nav_lat'      ,'f',('y','x')                        ); ncvar[:] = Lat;
        ncvar = ncOUT.createVariable('depthw'       ,'f',('depthw',)                      ); ncvar[:] = gdepw;
        ncvar = ncOUT.createVariable('time_counter' ,'d',('time_counter',)                ); ncvar    = 1.;
        ncvar = ncOUT.createVariable('vovecrtz'     ,'f',('time_counter','depthw','y','x')); ncvar[:] = D3*0.;
        ncvar = ncOUT.createVariable('votkeavt'     ,'f',('time_counter','depthw','y','x')); ncvar[:] = D3K;        

        ncOUT.close()
