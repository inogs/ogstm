import os,sys

import numpy as np

from mydtype import *

import scipy.io.netcdf as NC

import pickle

#import create_UVW as UVW
#import create_TSKQWHF  as TSKQWHF

def create_optics_nc(test):

    jpi=test['jpi'];
    jpj=test['jpj'];
    jpk=test['jpk'];
    time = 1
    maskfile=test['Dir'].decode() + '/meshmask.nc'

    M=NC.netcdf_file(maskfile,"r")

    Lon     =  M.variables['glamt'].data[0,0,:,:].copy()
    Lat     =  M.variables['gphit'].data[0,0,:,:].copy()
    gdept   =  M.variables['gdept'].data[0,:,0,0].copy()
    gdepw   =  M.variables['gdepw'].data[0,:,0,0].copy()
    e1v     =  M.variables['e1v'  ].data[0,0,:,:].copy()
    e2u     =  M.variables['e2u'  ].data[0,0,:,:].copy()

    e3v0    =  M.variables['e3v_0'  ].data[0,:,:,:].copy()
    e3u0    =  M.variables['e3u_0'  ].data[0,:,:,:].copy()
    e3t0    =  M.variables['e3t_0'  ].data[0,:,:,:].copy()
    e3w0    =  M.variables['e3w_0'  ].data[0,:,:,:].copy()
    
    
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



    OPTICS_DATE=[]

    filein=open('OPTICS/forcing_date_atm')
    for var in filein:
        OPTICS_DATE.append(var[:-1])

    filein.close()

    os.system("mkdir -p " + test['Dir'].decode() + '/OPTICS')

    for date in OPTICS_DATE:
        # Create T file
#       TSKQWHF.create_TSKQWHF(test,date,D3T,D3S,D3K,D2Q,D2W,D2H,D2F)
        yyyy=date[0:4]
        mm=date[4:6]
        outfile = test['Dir'].decode() + '/OPTICS/' + yyyy + '/' + mm + '/atm.' + date + '.nc'
        os.system("mkdir -p " + test['Dir'].decode() + '/OPTICS/' + yyyy + '/' + mm)
        ncOUT   = NC.netcdf_file(outfile,'w')

        ncOUT.createDimension('lon'           ,jpi);
        ncOUT.createDimension('lat'           ,jpj);

        ncvar = ncOUT.createVariable('lon'      ,'f',('lon',)                            ); ncvar[:] = Lon[0,:];
        ncvar = ncOUT.createVariable('lat'      ,'f',('lat',)                            ); ncvar[:] = Lat[:,0];
        ncvar = ncOUT.createVariable('sp'      ,'f',('lat','lon')                         ); ncvar[:] = 100000;
        ncvar = ncOUT.createVariable('msl'     ,'f',('lat','lon')                         ); ncvar[:] = 100000;
        ncvar = ncOUT.createVariable('u10'     ,'f',('lat','lon')          ); ncvar[:] = 5.0/np.sqrt(2.);
        ncvar = ncOUT.createVariable('v10'     ,'f',('lat','lon')          ); ncvar[:] = 5.0/np.sqrt(2.);
#       ncvar = ncOUT.createVariable('u10'     ,'f',('lat','lon')          ); ncvar[:] = D2W/np.sqrt(2.);
#       ncvar = ncOUT.createVariable('v10'     ,'f',('lat','lon')          ); ncvar[:] = D2W/np.sqrt(2.);
        ncvar = ncOUT.createVariable('tclw'     ,'f',('lat','lon')          ); ncvar[:] = 0.001;
        ncvar = ncOUT.createVariable('tco3'     ,'f',('lat','lon')          ); ncvar[:] = 0.005;
        ncvar = ncOUT.createVariable('t2m'     ,'f',('lat','lon')          ); ncvar[:] = 293.;
        ncvar = ncOUT.createVariable('d2m'     ,'f',('lat','lon')          ); ncvar[:] = 293.;
        ncvar = ncOUT.createVariable('tcc'     ,'f',('lat','lon')          ); ncvar[:] = 0.2;
#     	float lon(lon) ;
#	lon:units = "degrees" ;
#	lon:long_name = "longitude" ;
#float lat(lat) ;
#	lat:units = "degrees" ;
#	lat:long_name = "latitude" ;
#float sp(lat, lon) ;
#	sp:long_name = "surface pressure" ;
#	sp:units = "Pa" ;
#	sp:orig = "ERA5" ;
#float msl(lat, lon) ;
#	msl:long_name = "mean sea level pressure" ;
#	msl:units = "Pa" ;
#	msl:orig = "ERA5" ;
#float u10(lat, lon) ;
#	u10:long_name = "zonal wind velocity" ;
#	u10:units = "m s*-1" ;
#	u10:orig = "ERA5" ;
#float v10(lat, lon) ;
#	v10:long_name = "meridional wind velocity" ;
#	v10:units = "m s*-1" ;
#	v10:orig = "ERA5" ;
#float tclw(lat, lon) ;
#	tclw:long_name = "Total column cloud liquid water" ;
#	tclw:units = "kg m**-2" ;
#	tclw:orig = "ERA5" ;
#float tco3(lat, lon) ;
#	tco3:long_name = "Total column ozone" ;
#	tco3:units = "kg m**-2" ;
#	tco3:orig = "ERA5" ;
#float t2m(lat, lon) ;
#	t2m:long_name = "2 metre temperature" ;
#	t2m:units = "K" ;
#	t2m:orig = "ERA5" ;
#float d2m(lat, lon) ;
#	d2m:long_name = "2 metre dewpoint temperature" ;
#	d2m:units = "K" ;
#	d2m:orig = "ERA5" ;
#float tcc(lat, lon) ;
#	tcc:long_name = "Total cloud cover" ;
#	tcc:units = "[-]" ;
#	tcc:orig = "ERA5" ;
        ncOUT.close()

        # Create U file
    filein=open('OPTICS/forcing_date_atmclim')
    for var in filein:
        OPTICS_DATE.append(var[:-1])

    filein.close()

    os.system("mkdir -p " + test['Dir'].decode() + '/OPTICS')

    for date in OPTICS_DATE:

        outfile = test['Dir'].decode() + '/OPTICS/climatm.' + date + '.nc'
        ncOUT   = NC.netcdf_file(outfile,'w')


        ncOUT.createDimension('lon'           ,jpi);
        ncOUT.createDimension('lat'           ,jpj);

        ncvar = ncOUT.createVariable('lon'      ,'f',('lon',)                            ); ncvar[:] = Lon[0,:];
        ncvar = ncOUT.createVariable('lat'      ,'f',('lat',)                            ); ncvar[:] = Lat[:,0];

        ncvar = ncOUT.createVariable('cdrem'     ,'f',('lat','lon')          ); ncvar[:] = 12.;
        ncvar = ncOUT.createVariable('cldtcm'    ,'f',('lat','lon')          ); ncvar[:] = 15.;

        ncOUT.close()
#float cdrem(lat, lon) ;
#	cdrem:long_name = "cloud droplet effective radius" ;
#	cdrem:units = "[um]" ;
#	cdrem:orig = "MODCLD" ;
#float cldtcm(lat, lon) ;
#	cldtcm:long_name = "TODO" ;
#	cldtcm:units = "TODO" ;
#	cldtcm:orig = "MODCLD" ;

        # Create V file
    filein=open('OPTICS/forcing_date_aero')
    for var in filein:
        OPTICS_DATE.append(var[:-1])

    filein.close()

    os.system("mkdir -p " + test['Dir'].decode() + '/OPTICS')

    for date in OPTICS_DATE:
        outfile = test['Dir'].decode() + '/OPTICS/aero.' + date + '.nc'
        ncOUT   = NC.netcdf_file(outfile,'w')

        ncOUT.createDimension('lon'           ,jpi);
        ncOUT.createDimension('lat'           ,jpj);
        ncOUT.createDimension('depth'         ,jpk);

        ncvar = ncOUT.createVariable('lon'      ,'f',('lon',)                            ); ncvar[:] = Lon[0,:];
        ncvar = ncOUT.createVariable('lat'      ,'f',('lat',)                            ); ncvar[:] = Lat[:,0];

        ncvar = ncOUT.createVariable('taua'     ,'f',('depth','lat','lon')          ); ncvar[:] = 0.2;
        ncvar = ncOUT.createVariable('asymp'    ,'f',('depth','lat','lon')          ); ncvar[:] = 0.7;
        ncvar = ncOUT.createVariable('ssalb'    ,'f',('depth','lat','lon')          ); ncvar[:] = 0.98;

#float taua(depth, lat, lon) ;
#	taua:_FillValue = 1.e+20f ;
#	taua:long_name = "aerosol optical thickness" ;
#	taua:units = "[-]" ;
#	taua:orig = "MODIS_AEROSOL" ;
#float asymp(depth, lat, lon) ;
#	asymp:_FillValue = 1.e+20f ;
#	asymp:long_name = "aerosol asymmetry parameter" ;
#	asymp:units = "[-]" ;
#	asymp:orig = "MODIS_AEROSOL" ;
#float ssalb(depth, lat, lon) ;
#	ssalb:_FillValue = 1.e+20f ;
#	ssalb:long_name = "aerosol single scattering albedo" ;
#	ssalb:units = "[-]" ;
#	ssalb:orig = "MODIS_AEROSOL" ;

        ncOUT.close()

