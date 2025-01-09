import os,sys

import numpy as np

from mydtype import *

import scipy.io.netcdf as NC

import pickle

def create_bc_nc(test):

# Create boundaries namelist
    filename = test['Dir'].decode() + '/boundaries.nml'
    f01 = open(filename,'w')
    f01.write("3")
    f01.write("\n")
    f01.write("\"riv, RIV, riv.nml, files_namelist_riv.dat, T, F\"")
    f01.write("\"riv, SPO, gib.nml, files_namelist_gib.dat, T, F\"")
    f01.write("\"dar, OPE, dar.nml, files_namelist_dar.dat, T, F\"")
    f01.close()
# Create riv.nml namelist
    filename = test['Dir'].decode() + '/riv.nml'
    f01 = open(filename,'w')
    f01.write("&VARS_DIMENSION")
    f01.write("\n")
    f01.write("\n")
    f01.write("    n_vars = 6")
    f01.write("\n")
    f01.write("\n")
    f01.write("/")
    f01.write("\n")
    f01.write("\n")
    f01.write("&CORE")
    f01.write("\n")
    f01.write("\n")
    f01.write("    vars(1) = \"N1p\"")
    f01.write("\n")
    f01.write("    vars(2) = \"N3n\"")
    f01.write("\n")
    f01.write("    vars(3) = \"N5s\"")
    f01.write("\n")
    f01.write("    vars(4) = \"O3c\"")
    f01.write("\n")
    f01.write("    vars(5) = \"O3h\"")
    f01.write("\n")
    f01.write("    vars(6) = \"O2o\"")
    f01.write("\n")
    f01.write("\n")
    f01.write("/")
    f01.write("\n")
    f01.close()
# Create files_namelist_riv.dat namelist
    filename = test['Dir'].decode() + '/files_namelist_riv.dat'
    f01 = open(filename,'w')
    f01.write("12")
    f01.write("\n")
    f01.write("\"BC/TIN_yyyy0115-00:00:00.nc\"")
    f01.write("\n")
    f01.write("\"BC/TIN_yyyy0215-00:00:00.nc\"")
    f01.write("\n")
    f01.write("\"BC/TIN_yyyy0315-00:00:00.nc\"")
    f01.write("\n")
    f01.write("\"BC/TIN_yyyy0415-00:00:00.nc\"")
    f01.write("\n")
    f01.write("\"BC/TIN_yyyy0515-00:00:00.nc\"")
    f01.write("\n")
    f01.write("\"BC/TIN_yyyy0615-00:00:00.nc\"")
    f01.write("\n")
    f01.write("\"BC/TIN_yyyy0715-00:00:00.nc\"")
    f01.write("\n")
    f01.write("\"BC/TIN_yyyy0815-00:00:00.nc\"")
    f01.write("\n")
    f01.write("\"BC/TIN_yyyy0915-00:00:00.nc\"")
    f01.write("\n")
    f01.write("\"BC/TIN_yyyy1015-00:00:00.nc\"")
    f01.write("\n")
    f01.write("\"BC/TIN_yyyy1115-00:00:00.nc\"")
    f01.write("\n")
    f01.write("\"BC/TIN_yyyy1215-00:00:00.nc\"")
    f01.close()

    filename = test['Dir'].decode() + '/files_namelist_gib.dat'
    f01 = open(filename,'w')
    f01.write("4")
    f01.write("\n")
    f01.write("\"BC/GIB_yyyy0215-00:00:00.nc\"")
    f01.write("\n")
    f01.write("\"BC/GIB_yyyy0515-00:00:00.nc\"")
    f01.write("\n")
    f01.write("\"BC/GIB_yyyy0815-00:00:00.nc\"")
    f01.write("\n")
    f01.write("\"BC/GIB_yyyy1115-00:00:00.nc\"")
    f01.write("\n")
    f01.close()

    filename = test['Dir'].decode() + '/files_namelist_dar.dat'
    f01 = open(filename,'w')
    f01.write("1\n")
    f01.write("\"BC/OPE_yyyy0630-00:00:00.nc\"")
    f01.write("\n")
    f01.close()

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
    tmask   =  M.variables['tmask'].data[0,:,:,:].astype(bool).copy()
    M.close()

    D3=np.ones((1,jpk,jpj,jpi),float)
    D2=np.ones((1,jpj,jpi),float)

# Creating bounmask.nc

    index=np.zeros((jpk,jpj,jpi),int)

    waterpoints = 0;

    for jk in range(jpk):
        for jj in range(jpj):
            for ji in range(jpi):
                if tmask[jk,jj,ji]:
                    waterpoints += 1
                    index[jk,jj,ji] = waterpoints

    index_inv=np.zeros((waterpoints,3),int)

    waterpoints = 0;

    for jk in range(jpk):
        for jj in range(jpj):
            for ji in range(jpi):
                if tmask[jk,jj,ji]:
                    index_inv[waterpoints,0] = jk + 1 # Fortran
                    index_inv[waterpoints,1] = jj + 1
                    index_inv[waterpoints,2] = ji + 1
                    waterpoints += 1

    outfile = test['Dir'].decode() + '/bounmask.nc'
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

    os.system("mkdir -p " + test['Dir'].decode() + '/BC/')
# Atmosphere NUT
    ATM_DATE=[]

    filein=open('KB/atm_date')

    for var in filein:
        ATM_DATE.append(var[:-1])

    filein.close()


    atm_idxt = 0;

    for jj in range(jpj):
        for ji in range(jpi):
                if tmask[0,jj,ji]:
                    atm_idxt += 1

    atm_index=np.ones((atm_idxt),float)

    atm_idxt = 0;

    for jj in range(jpj):
        for ji in range(jpi):
                if tmask[0,jj,ji]:
                    atm_index[atm_idxt] = index[0,jj,ji]
                    atm_idxt            += 1

    for date in ATM_DATE:
        # Create ATM file
        outfile = test['Dir'].decode() + '/BC/ATM_' + date + '.nc'
        ncOUT   = NC.netcdf_file(outfile,'w')


        ncOUT.createDimension('lon'    ,jpi);
        ncOUT.createDimension('lat'    ,jpj);        
        ncvar = ncOUT.createVariable('atm_N1p'      ,'f',('lat','lon')                   ); ncvar[:] = 3.75866672509673e-09;
        ncvar = ncOUT.createVariable('atm_N3n'      ,'f',('lat','lon')                   ); ncvar[:] = 2.24183651189621e-07;

        ncOUT.close()

# Atmosphere CO2

    CO2_DATE=[]

    filein=open('KB/co2_date')

    for var in filein:
        CO2_DATE.append(var[:-1])

    filein.close()

    D2=np.ones((jpj,jpi),float)


    for date in CO2_DATE:
        # Create CO2 file
        outfile = test['Dir'].decode() + '/BC/CO2_' + date + '.nc'
        ncOUT   = NC.netcdf_file(outfile,'w')

        ncOUT.createDimension('lon',jpi);
        ncOUT.createDimension('lat',jpj);

        ncvar = ncOUT.createVariable('CO2'     ,'f',('lat','lon')                   ); ncvar[:] = 390.;

        ncOUT.close()

# TIN
    TIN_DATE=[]

    filein=open('KB/tin_date')

    for var in filein:
        TIN_DATE.append(var[:-1])

    filein.close()



    tin_idxt = 0;

    for wp in range(waterpoints):
        ji = index_inv[wp,2] -1 ; # index_inv is fortran style
        jj = index_inv[wp,1] -1 ;
        jk = index_inv[wp,0] -1 ;
        if ( (ji==1) | (jj==1) ) | ( (ji==jpi-2) | (jj==jpj-2) ):
            if (jk ==0 ):
                tin_idxt += 1

    riv_index=np.ones((tin_idxt),int)

    tin_idxt = 0;

    for wp in range(waterpoints):
        ji = index_inv[wp,2] -1; # index_inv is fortran style
        jj = index_inv[wp,1] -1;
        jk = index_inv[wp,0] -1;
        if ( (ji==1) | (jj==1) ) | ( (ji==jpi-2) | (jj==jpj-2) ):
            if (jk ==0 ):
                riv_index[tin_idxt]=index[jk,jj,ji]
                tin_idxt += 1

    riv_N1p=np.zeros((jpj,jpi),dtype=float)-1.0; riv_N1p[0,0]=0.35*10**(-5)
    riv_N3n=np.zeros((jpj,jpi),dtype=float)-1.0; riv_N3n[0,0]=0.2*10**(-3)
    riv_N5s=np.zeros((jpj,jpi),dtype=float)-1.0; riv_N5s[0,0]=1.0*10**(-4)
    riv_O3c=np.zeros((jpj,jpi),dtype=float)-1.0; riv_O3c[0,0]=0.35
    riv_O3h=np.zeros((jpj,jpi),dtype=float)-1.0; riv_O3h[0,0]=0.01
    riv_O2o=np.zeros((jpj,jpi),dtype=float)-1.0; riv_O2o[0,0]=0.002

    for date in TIN_DATE:
        # Create RIV file
        outfile = test['Dir'].decode() + '/BC/TIN_' + date + '.nc'
        ncOUT   = NC.netcdf_file(outfile,'w')

        ncOUT.createDimension('lon'    ,jpi);
        ncOUT.createDimension('lat'    ,jpj);

        ncvar = ncOUT.createVariable('riv_N1p'      ,'d',('lat','lon') ); ncvar[:] = riv_N1p;
        ncvar.missing_value = 1e20
        ncvar = ncOUT.createVariable('riv_N3n'      ,'d',('lat','lon') ); ncvar[:] = riv_N3n;
        ncvar.missing_value = 1e20
        ncvar = ncOUT.createVariable('riv_N5s'      ,'d',('lat','lon') ); ncvar[:] = riv_N5s;
        ncvar.missing_value = 1e20
        ncvar = ncOUT.createVariable('riv_O3c'      ,'d',('lat','lon') ); ncvar[:] = riv_O3c;
        ncvar.missing_value = 1e20
        ncvar = ncOUT.createVariable('riv_O3h'      ,'d',('lat','lon') ); ncvar[:] = riv_O3h;
        ncvar.missing_value = 1e20
        ncvar = ncOUT.createVariable('riv_O2o'      ,'d',('lat','lon') ); ncvar[:] = riv_O2o;
        ncvar.missing_value = 1e20
        ncOUT.close()

# GIB
    GIB_DATE=[]

    filein=open('KB/gib_date')

    for var in filein:
        GIB_DATE.append(var[:-1])

    filein.close()



    gib_idxt = 0;

    for wp in range(waterpoints):
        ji = index_inv[wp,2] -1 ; # index_inv is fortran style
        jj = index_inv[wp,1] -1 ;
        if ( (ji==1) | (jj==1) ) | ( (ji==jpi-2) | (jj==jpj-2) ):
            gib_idxt += 1

    gib_index=np.ones((gib_idxt),int)

    gib_idxt = 0;

    for wp in range(waterpoints):
        ji = index_inv[wp,2] -1 ; # index_inv is fortran style
        jj = index_inv[wp,1] -1 ;
        jk = index_inv[wp,0] -1 ; 
        if ( (ji==1) | (jj==1) ) | ( (ji==jpi-2) | (jj==jpj-2) ):
            gib_index[gib_idxt]=index[jk,jj,ji]
            gib_idxt += 1


    for date in GIB_DATE:
        # Create GIB file
        outfile = test['Dir'].decode() + '/BC/GIB_' + date + '.nc'
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

        BOUNDARY_CONCENTRATION={}
        BOUNDARY_CONCENTRATION['N5s'] = 2.0
        BOUNDARY_CONCENTRATION['N1p'] = 0.065 # mmol/m3
        BOUNDARY_CONCENTRATION['N3n']=   1.3 # mol/m3
        BOUNDARY_CONCENTRATION['O3c']= 28700 # mg/m3
        BOUNDARY_CONCENTRATION['O3h']=  2800 # mmol/m3
        OPEN_BOUNDARY=np.zeros((jpk,jpj,jpi),np.float32)
        I = jpi-1
        j_min=1
        j_max=int(jpj/2)
        for k in range(jpk):
            OPEN_BOUNDARY[k,j_min:j_max,I] = 0.1
        outfile = test['Dir'].decode() + '/BC/OPE_yyyy0630-00:00:00.nc'
        ncOUT   = NC.netcdf_file(outfile,'w')
        ncOUT.createDimension('lon'    ,jpi);
        ncOUT.createDimension('lat'    ,jpj);
        ncOUT.createDimension('z'     ,jpk );
        for var in OPEN_BOUNDARY.keys():
            for k in range(jpk):
                OPEN_BOUNDARY[k,j_min:j_max,I] = BOUNDARY_CONCENTRATION[var]
            ncvar = ncOUT.createVariable(var,f,('z','lat','lon')  )
            ncvar[:] = OPEN_BOUNDARY

        ncOUT.close()

