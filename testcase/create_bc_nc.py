import os,sys

import numpy as np

from mydtype import *

import scipy.io.netcdf as NC

import pickle

def create_bc_nc(test):

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
    tmask   =  M.variables['tmask'].data[0,:,:,:].astype(bool).copy()
    M.close()

    D3=np.ones((1,jpk,jpj,jpi),np.float)    
    D2=np.ones((1,jpj,jpi),np.float)

# Creating bounmask.nc

    index=np.zeros((jpk,jpj,jpi),np.int)

    waterpoints = 0;

    for jk in range(jpk):
        for jj in range(jpj):
            for ji in range(jpi):
                if tmask[jk,jj,ji]:
                    waterpoints += 1
                    index[jk,jj,ji] = waterpoints

    index_inv=np.zeros((waterpoints,3),np.int)

    waterpoints = 0;

    for jk in range(jpk):
        for jj in range(jpj):
            for ji in range(jpi):
                if tmask[jk,jj,ji]:
                    index_inv[waterpoints,0] = jk + 1 # Fortran  
                    index_inv[waterpoints,1] = jj + 1
                    index_inv[waterpoints,2] = ji + 1
                    waterpoints += 1

    outfile = test['Dir'] + '/bounmask.nc'
    ncOUT   = NC.netcdf_file(outfile,'w')
    ncOUT.createDimension('x'          ,jpi );
    ncOUT.createDimension('y'          ,jpj ); 
    ncOUT.createDimension('z'          ,jpk );
    ncOUT.createDimension('time'       ,time);
    ncOUT.createDimension('waterpoints',waterpoints   ); # Points to considered are on the vertex
    ncOUT.createDimension('dim3'       ,3   );

    ncvar = ncOUT.createVariable('nav_lon'      ,'f',('y','x')             ); ncvar[:] = Lon;
    ncvar = ncOUT.createVariable('nav_lat'      ,'f',('y','x')             ); ncvar[:] = Lat;
    ncvar = ncOUT.createVariable('nav_lev'      ,'f',('z',)                ); ncvar[:] = gdept;

    ncvar = ncOUT.createVariable('reN1p'        ,'d',('time','z','y','x')  ); ncvar[:] = D3*0.;   
    ncvar = ncOUT.createVariable('reN3n'        ,'d',('time','z','y','x')  ); ncvar[:] = D3*0.;
    ncvar = ncOUT.createVariable('reO2o'        ,'d',('time','z','y','x')  ); ncvar[:] = D3*0.;
    ncvar = ncOUT.createVariable('reN5s'        ,'d',('time','z','y','x')  ); ncvar[:] = D3*0.;
    ncvar = ncOUT.createVariable('reO3c'        ,'d',('time','z','y','x')  ); ncvar[:] = D3*0.; 
    ncvar = ncOUT.createVariable('reO3h'        ,'d',('time','z','y','x')  ); ncvar[:] = D3*0.;
    ncvar = ncOUT.createVariable('reN6r'        ,'d',('time','z','y','x')  ); ncvar[:] = D3*0.;

    ncvar = ncOUT.createVariable('index'        ,'i',('time','z','y','x')  ); ncvar[:] = index;
    ncvar = ncOUT.createVariable('index_inv'    ,'i',('waterpoints','dim3')); ncvar[:] = index_inv;

    ncOUT.close()

    os.system("mkdir -p " + test['Dir'] + '/BC/')
# Atmosphere NUT
    ATM_DATE=[]

    filein=file('KB/atm_date')

    for var in filein:
        ATM_DATE.append(var[:-1])
    
    filein.close()    


    atm_idxt = 0;

    for jj in range(jpj):
        for ji in range(jpi):
                if tmask[0,jj,ji]:
                    atm_idxt += 1

    atm_index=np.ones((atm_idxt),np.float)

    atm_idxt = 0;

    for jj in range(jpj):
        for ji in range(jpi):
                if tmask[0,jj,ji]:
                    atm_index[atm_idxt] = index[0,jj,ji]
                    atm_idxt            += 1

    for date in ATM_DATE:
        # Create ATM file
        outfile = test['Dir'] + '/BC/ATM_' + date + '.nc'
        ncOUT   = NC.netcdf_file(outfile,'w')

        ncOUT.createDimension('atm_idxt'    ,atm_idxt);
        
        ncvar = ncOUT.createVariable('atm_idxt'     ,'i',('atm_idxt',)                   ); ncvar[:] = atm_index;
        ncvar = ncOUT.createVariable('atm_N1p'      ,'d',('atm_idxt',)                   ); ncvar[:] = 3.75866672509673e-09;
        ncvar = ncOUT.createVariable('atm_N3n'      ,'d',('atm_idxt',)                   ); ncvar[:] = 2.24183651189621e-07;
        ncOUT.close()

# Atmosphere CO2
 
    CO2_DATE=[]

    filein=file('KB/co2_date')

    for var in filein:
        CO2_DATE.append(var[:-1])
    
    filein.close()    

    D2=np.ones((jpj,jpi),np.float)


    for date in CO2_DATE:
        # Create CO2 file
        outfile = test['Dir'] + '/BC/CO2_' + date + '.nc'
        ncOUT   = NC.netcdf_file(outfile,'w')

        ncOUT.createDimension('lon',jpi);
        ncOUT.createDimension('lat',jpj);
        
        ncvar = ncOUT.createVariable('CO2'     ,'f',('lat','lon')                   ); ncvar[:] = 390.;

        ncOUT.close()
 
# TIN 
    TIN_DATE=[]

    filein=file('KB/tin_date')

    for var in filein:
        TIN_DATE.append(var[:-1])
    
    filein.close()    

    

    tin_idxt = 0;

    for wp in range(waterpoints):
        ji = index_inv[wp,2] -1 ; # index_inv is fortran style
        jj = index_inv[wp,1] -1 ;
        jk = index_inv[wp,0] -1 ;
        if ( (ji==1) | (jj==1) ) | ( (ji==jpi-2) | (jj==jpj-2) ):
		if (jk < 2):
            		tin_idxt += 1

    riv_index=np.ones((tin_idxt),np.int)

    tin_idxt = 0;

    for wp in range(waterpoints):
        ji = index_inv[wp,2] -1; # index_inv is fortran style
        jj = index_inv[wp,1] -1;
        jk = index_inv[wp,0] -1;
        if ( (ji==1) | (jj==1) ) | ( (ji==jpi-2) | (jj==jpj-2) ):
		if (jk < 2):
            		riv_index[tin_idxt]=index[jk,jj,ji]
            		tin_idxt += 1

    riv_N1p=np.zeros(tin_idxt,dtype=float); riv_N1p[0]=0.35*10**(-5)
    riv_N3n=np.zeros(tin_idxt,dtype=float); riv_N3n[0]=0.2*10**(-3)
    riv_N5s=np.zeros(tin_idxt,dtype=float); riv_N5s[0]=1.0*10**(-4)
    riv_O3c=np.zeros(tin_idxt,dtype=float); riv_O3c[0]=0.35
    riv_O3h=np.zeros(tin_idxt,dtype=float); riv_O3h[0]=0.01

    for date in TIN_DATE:
        # Create RIV file
        outfile = test['Dir'] + '/BC/TIN_' + date + '.nc'
        ncOUT   = NC.netcdf_file(outfile,'w')

        ncOUT.createDimension('riv_idxt'    ,tin_idxt);
        
        ncvar = ncOUT.createVariable('riv_idxt'     ,'i',('riv_idxt',) ); ncvar[:] = riv_index;
        ncvar = ncOUT.createVariable('riv_N1p'      ,'d',('riv_idxt',) ); ncvar[:] = riv_N1p;
        ncvar = ncOUT.createVariable('riv_N3n'      ,'d',('riv_idxt',) ); ncvar[:] = riv_N3n;
        ncvar = ncOUT.createVariable('riv_N5s'      ,'d',('riv_idxt',) ); ncvar[:] = riv_N5s;
        ncvar = ncOUT.createVariable('riv_O3c'      ,'d',('riv_idxt',) ); ncvar[:] = riv_O3c;
        ncvar = ncOUT.createVariable('riv_O3h'      ,'d',('riv_idxt',) ); ncvar[:] = riv_O3h;

        ncOUT.close()

# GIB 
    GIB_DATE=[]

    filein=file('KB/gib_date')

    for var in filein:
        GIB_DATE.append(var[:-1])
    
    filein.close()    

    

    gib_idxt = 0;

    for wp in range(waterpoints):
        ji = index_inv[wp,2] -1 ; # index_inv is fortran style
        jj = index_inv[wp,1] -1 ;
        if ( (ji==1) | (jj==1) ) | ( (ji==jpi-2) | (jj==jpj-2) ):
            gib_idxt += 1

    gib_index=np.ones((gib_idxt),np.int)

    gib_idxt = 0;

    for wp in range(waterpoints):
        ji = index_inv[wp,2] -1 ; # index_inv is fortran style
        jj = index_inv[wp,1] -1 ;
        if ( (ji==1) | (jj==1) ) | ( (ji==jpi-2) | (jj==jpj-2) ):
            gib_index[gib_idxt]=index[0,jj,ji]
            gib_idxt += 1


    for date in GIB_DATE:
        # Create GIB file
        outfile = test['Dir'] + '/BC/GIB_' + date + '.nc'
        ncOUT   = NC.netcdf_file(outfile,'w')

        ncOUT.createDimension('gib_idxt_N1p'    ,gib_idxt);
        ncOUT.createDimension('gib_idxt_N3n'    ,gib_idxt);
        ncOUT.createDimension('gib_idxt_O2o'    ,gib_idxt);
        ncOUT.createDimension('gib_idxt_N5s'    ,gib_idxt);
        ncOUT.createDimension('gib_idxt_O3c'    ,gib_idxt);
        ncOUT.createDimension('gib_idxt_O3h'    ,gib_idxt);
        ncOUT.createDimension('gib_idxt_N6r'    ,gib_idxt);


       
        ncvar = ncOUT.createVariable('gib_idxt_N1p','i',('gib_idxt_N1p',) ); ncvar[:] = gib_index;
        ncvar = ncOUT.createVariable('gib_idxt_N3n','i',('gib_idxt_N3n',) ); ncvar[:] = gib_index;
        ncvar = ncOUT.createVariable('gib_idxt_O2o','i',('gib_idxt_O2o',) ); ncvar[:] = gib_index;
        ncvar = ncOUT.createVariable('gib_idxt_N5s','i',('gib_idxt_N5s',) ); ncvar[:] = gib_index;
        ncvar = ncOUT.createVariable('gib_idxt_O3c','i',('gib_idxt_O3c',) ); ncvar[:] = gib_index;
        ncvar = ncOUT.createVariable('gib_idxt_O3h','i',('gib_idxt_O3h',) ); ncvar[:] = gib_index;
        ncvar = ncOUT.createVariable('gib_idxt_N6r','i',('gib_idxt_N6r',) ); ncvar[:] = gib_index;

        ncvar = ncOUT.createVariable('gib_N1p'     ,'d',('gib_idxt_N1p',) ); ncvar[:] = 0.;
        ncvar = ncOUT.createVariable('gib_N3n'     ,'d',('gib_idxt_N3n',) ); ncvar[:] = 0.;
        ncvar = ncOUT.createVariable('gib_O2o'     ,'d',('gib_idxt_O2o',) ); ncvar[:] = 0.;
        ncvar = ncOUT.createVariable('gib_N5s'     ,'d',('gib_idxt_N5s',) ); ncvar[:] = 0.;
        ncvar = ncOUT.createVariable('gib_O3c'     ,'d',('gib_idxt_O3c',) ); ncvar[:] = 0.;
        ncvar = ncOUT.createVariable('gib_O3h'     ,'d',('gib_idxt_O3h',) ); ncvar[:] = 0.;
        ncvar = ncOUT.createVariable('gib_N6r'     ,'d',('gib_idxt_N6r',) ); ncvar[:] = 0.;

        ncOUT.close()
