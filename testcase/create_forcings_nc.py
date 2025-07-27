import os,sys

import numpy as np

from mydtype import *

import  netCDF4 as NC

import pickle

import create_UVW as UVW
import create_TSKQWHF  as TSKQWHF

def create_forcings_nc(test):

    jpi=test['jpi'];
    jpj=test['jpj'];
    jpk=test['jpk'];
    time = 1
    maskfile=test['Dir'].decode() + '/meshmask.nc'

    M=NC.Dataset(maskfile,"r")

    Lon     =  M.variables['glamt'][0,0,:,:].copy()
    Lat     =  M.variables['gphit'][0,0,:,:].copy()
    gdept   =  M.variables['gdept'][0,:,0,0].copy()
    gdepw   =  M.variables['gdepw'][0,:,0,0].copy()
    e1v     =  M.variables['e1v'  ][0,0,:,:].copy()
    e2u     =  M.variables['e2u'  ][0,0,:,:].copy()

    e3v0    =  M.variables['e3v_0'  ][0,:,:,:].copy()
    e3u0    =  M.variables['e3u_0'  ][0,:,:,:].copy()
    e3t0    =  M.variables['e3t_0'  ][0,:,:,:].copy()
    e3w0    =  M.variables['e3w_0'  ][0,:,:,:].copy()
    
    
    M.close()
    e3t = e3t0
    e1v_aux = np.tile(e1v,[jpk,1,1])
    e2u_aux = np.tile(e2u,[jpk,1,1])

    Av=np.ones((jpk,jpj,jpi),np.float64)
    Au=np.ones((jpk,jpj,jpi),np.float64)

    for i in range(jpk):
        Av[i,:,:] = e1v_aux[i,:,:] * e3t[i]
        Au[i,:,:] = e2u_aux[i,:,:] * e3t[i]


    D3T=np.zeros((1,jpk,jpj,jpi),np.float64)
    D3S=np.zeros((1,jpk,jpj,jpi),np.float64)
    D3K=np.zeros((1,jpk,jpj,jpi),np.float64)

    D2Q=np.zeros((1,jpj,jpi),np.float64)
    D2W=np.zeros((1,jpj,jpi),np.float64)
    D2H=np.zeros((1,jpj,jpi),np.float64)
    D2F=np.zeros((1,jpj,jpi),np.float64)

    D3U=np.zeros((1,jpk,jpj,jpi),np.float64)
    D3V=np.zeros((1,jpk,jpj,jpi),np.float64)
    D3W=np.zeros((1,jpk,jpj,jpi),np.float64)

    UVW.create_UVW(test,D3U,D3V,D3W,Av,Au)

    D3=np.ones((1,jpk,jpj,jpi),np.float64)
    D2=np.ones((1,jpj,jpi),np.float64)
    SSH=np.zeros((1,jpj,jpi),np.float64)-0.6


    FORCING_DATE=[]

    filein=open('KB/forcing_date')
    for var in filein:
        FORCING_DATE.append(var[:-1])

    filein.close()

    os.system("mkdir -p " + test['Dir'].decode() + '/FORCINGS/')
    os.system("mkdir -p " + test['Dir'].decode() + '/FORCINGS/yyyy')

    for date in FORCING_DATE:
        yyyy=date[0:4]
        mm=date[4:6]
        # Create T file
        TSKQWHF.create_TSKQWHF(test,date,D3T,D3S,D3K,D2Q,D2W,D2H,D2F)
        outfile = test['Dir'].decode() + '/FORCINGS/' + yyyy + '/' + mm + '/T' + date + '.nc'
        os.system("mkdir -p " + test['Dir'].decode() + '/FORCINGS/yyyy/' + mm)
        ncOUT   = NC.Dataset(outfile,'w')

        ncOUT.createDimension('time_counter',None);
        ncOUT.createDimension('x'           ,jpi);
        ncOUT.createDimension('y'           ,jpj);
        ncOUT.createDimension('deptht'      ,jpk);

        ncvar = ncOUT.createVariable('nav_lon'      ,'f',('y','x')                        ); ncvar[:] = Lon;
        ncvar = ncOUT.createVariable('nav_lat'      ,'f',('y','x')                        ); ncvar[:] = Lat;
        ncvar = ncOUT.createVariable('deptht'       ,'f',('deptht',)                      ); ncvar[:] = gdept;
        ncvar = ncOUT.createVariable('time_counter' ,'d',('time_counter',)                ); ncvar    = 1.;
        ncvar = ncOUT.createVariable('vosaline'     ,'f',('time_counter','deptht','y','x')); ncvar[:] = D3S;
        ncvar = ncOUT.createVariable('votemper'     ,'f',('time_counter','deptht','y','x')); ncvar[:] = D3T;
        ncvar = ncOUT.createVariable('soshfldo'     ,'f',('time_counter','y','x')         ); ncvar[:] = D2Q;
        ncvar = ncOUT.createVariable('sowindsp'     ,'f',('time_counter','y','x')         ); ncvar[:] = D2W;
        ncvar = ncOUT.createVariable('sossheig'     ,'f',('time_counter','y','x')         ); ncvar[:] = SSH;
         
        ncOUT.close()

        # Create U file

        outfile = test['Dir'].decode() + '/FORCINGS/' + yyyy + '/' + mm + '/U' + date + '.nc'
        ncOUT   = NC.Dataset(outfile,'w')

        ncOUT.createDimension('time_counter',None);
        ncOUT.createDimension('x'           ,jpi);
        ncOUT.createDimension('y'           ,jpj);
        ncOUT.createDimension('depthu'      ,jpk);

        ncvar = ncOUT.createVariable('nav_lon'      ,'f',('y','x')                        ); ncvar[:] = Lon;
        ncvar = ncOUT.createVariable('nav_lat'      ,'f',('y','x')                        ); ncvar[:] = Lat;
        ncvar = ncOUT.createVariable('depthu'       ,'f',('depthu',)                      ); ncvar[:] = gdept;
        ncvar = ncOUT.createVariable('time_counter' ,'d',('time_counter',)                ); ncvar    = 1.;
        ncvar = ncOUT.createVariable('vozocrtx'     ,'f',('time_counter','depthu','y','x')); ncvar[:] = D3U;  
        ncvar = ncOUT.createVariable('sozotaux'     ,'f',('time_counter','y','x')); ncvar[:] = D2;  

        ncOUT.close()

        # Create V file

        outfile = test['Dir'].decode() + '/FORCINGS/' + yyyy + '/' + mm + '/V' + date + '.nc'
        ncOUT   = NC.Dataset(outfile,'w')

        ncOUT.createDimension('time_counter',None);
        ncOUT.createDimension('x'           ,jpi);
        ncOUT.createDimension('y'           ,jpj);
        ncOUT.createDimension('depthv'      ,jpk);

        ncvar = ncOUT.createVariable('nav_lon'      ,'f',('y','x')                        ); ncvar[:] = Lon;
        ncvar = ncOUT.createVariable('nav_lat'      ,'f',('y','x')                        ); ncvar[:] = Lat;
        ncvar = ncOUT.createVariable('depthv'       ,'f',('depthv',)                      ); ncvar[:] = gdept;
        ncvar = ncOUT.createVariable('time_counter'  ,'d',('time_counter',)               ); ncvar    = 1.;
        ncvar = ncOUT.createVariable('vomecrty'     ,'f',('time_counter','depthv','y','x')); ncvar[:] = D3V; 
        ncvar = ncOUT.createVariable('sometauy'     ,'f',('time_counter','y','x')); ncvar[:] = D2;  
        ncOUT.close()

        # Create W file

        outfile = test['Dir'].decode() + '/FORCINGS/' + yyyy + '/' + mm + '/W' + date + '.nc'
        ncOUT   = NC.Dataset(outfile,'w')

        ncOUT.createDimension('time_counter',None);
        ncOUT.createDimension('x'           ,jpi);
        ncOUT.createDimension('y'           ,jpj);
        ncOUT.createDimension('depthw'      ,jpk);

        ncvar = ncOUT.createVariable('nav_lon'      ,'f',('y','x')                        ); ncvar[:] = Lon;
        ncvar = ncOUT.createVariable('nav_lat'      ,'f',('y','x')                        ); ncvar[:] = Lat;
        ncvar = ncOUT.createVariable('depthw'       ,'f',('depthw',)                      ); ncvar[:] = gdepw;
        ncvar = ncOUT.createVariable('time_counter' ,'d',('time_counter',)                ); ncvar    = 1.;
        ncvar = ncOUT.createVariable('vovecrtz'     ,'f',('time_counter','depthw','y','x')); ncvar[:] = D3*0.;
        ncvar = ncOUT.createVariable('votkeavt'     ,'f',('time_counter','depthw','y','x')); ncvar[:] = D3K;               


        ncOUT.close()
