import os,sys
import numpy as np
from mydtype import test_conf, ext_data, file2stringlist
import scipy.io.netcdf as NC


def create_bc_nc(test):

# Create boundaries namelist
    filename = test['Dir'] + '/boundaries.nml'
    f01 = open(filename,'w')
    f01.write("3")
    f01.write("\n")
    f01.write("\"riv, RIV, riv.nml, files_namelist_riv.dat, T, F\"\n")
    f01.write("\"gib, SPO, gib.nml, files_namelist_gib.dat, T, T\"\n")
    f01.write("\"dar, OPE, dar.nml, files_namelist_dar.dat, T, F\"\n")
    f01.close()
# Create riv.nml namelist
    filename = test['Dir'] + '/riv.nml'
    f01 = open(filename,'w')
    f01.write("&VARS_DIMENSION\n\n")
    f01.write("    n_vars = 6\n\n")
    f01.write("/\n\n")
    f01.write("&CORE\n\n")
    f01.write("    vars(1) = \"N1p\"\n")
    f01.write("    vars(2) = \"N3n\"\n")
    f01.write("    vars(3) = \"N5s\"\n")
    f01.write("    vars(4) = \"O3c\"\n")
    f01.write("    vars(5) = \"O3h\"\n")
    f01.write("    vars(6) = \"O2o\"\n")
    f01.write("\n")
    f01.write("/")
    f01.write("\n")
    f01.close()

    TIN_DATE=file2stringlist('KB/tin_date')
    filename = test['Dir'] + '/files_namelist_riv.dat'
    f01 = open(filename,'w')
    f01.write(f"{len(TIN_DATE)}\n")
    for d in TIN_DATE:
        f01.write(f"\"BC/TIN_{d}.nc\"\n")
    f01.close()

    GIB_DATE=file2stringlist('KB/gib_date')
    filename = test['Dir'] + '/files_namelist_gib.dat'
    f01 = open(filename,'w')
    f01.write(f"{len(GIB_DATE)}\n")
    for d in GIB_DATE:
        f01.write(f"\"BC/GIB_{d}.nc\"\n")
    f01.close()

    filename = test['Dir'] + '/gib.nml'
    f01 = open(filename,'w')
    f01.write("&VARS_DIMENSION\n")
    f01.write("n_vars = 2\n")
    f01.write("/\n")
    f01.write("\n")
    f01.write("&CORE\n")
    f01.write("vars(1) = \"N1p\"\n")
    f01.write("vars(2) = \"N3n\"\n")
    f01.write("alpha = 4.0d0\n")
    f01.write("reduction_value_t = 1.0d-6\n")
    f01.write("length = -7.5d0\n")
    f01.write("/\n")
    f01.write("\n")
    f01.write("&NUDGING_VARS_DIMENSION\n")
    f01.write("n_vars = 2\n")
    f01.write("/\n")
    f01.write("&NUDGING_CORE\n")
    f01.write("data_file=\"bounmask.nc\"\n")
    f01.write("vars(1) = \"N1p\"\n")
    f01.write("vars(2) = \"N3n\"\n")
    f01.write("rst_corr(1) = 1.0d0\n")
    f01.write("rst_corr(2) = 1.0d0\n")
    f01.write("/\n")
    f01.close()
    
    
    filename = test['Dir'] + '/files_namelist_dar.dat'
    f01 = open(filename,'w')
    f01.write("1\n")
    f01.write("\"BC/OPE_yyyy0630-00:00:00.nc\"")
    f01.write("\n")
    f01.close()

    filename=test['Dir'] + '/dar.nml'
    f01 = open(filename,'w')
    f01.write("&VARS_DIMENSION\n")
    f01.write("n_vars = 5\n")
    f01.write("/\n")
    f01.write("\n")
    f01.write("&CORE\n")
    f01.write("vars(1) = \"N1p\"\n")
    f01.write("vars(2) = \"N3n\"\n")
    f01.write("vars(3) = \"N5s\"\n")
    f01.write("vars(4) = \"O3c\"\n")
    f01.write("vars(5) = \"O3h\"\n")
    f01.write("geometry = 1\n")
    f01.write("damping_coeff = 600.0d0\n")
    f01.write("/\n")
    f01.close()


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

    D3=np.zeros((1,jpk,jpj,jpi),float)
    D2=np.ones((1,jpj,jpi),float)
    first_quarter=int(jpi/4)
    D3[0,:,:,0:first_quarter]=1.e-05 # resto value

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

    ncvar = ncOUT.createVariable('reN1p'        ,'d',('time','z','y','x')  ); ncvar[:] = D3
    ncvar = ncOUT.createVariable('reN3n'        ,'d',('time','z','y','x')  ); ncvar[:] = D3
    ncvar = ncOUT.createVariable('reO2o'        ,'d',('time','z','y','x')  ); ncvar[:] = D3
    ncvar = ncOUT.createVariable('reN5s'        ,'d',('time','z','y','x')  ); ncvar[:] = D3
    ncvar = ncOUT.createVariable('reO3c'        ,'d',('time','z','y','x')  ); ncvar[:] = D3
    ncvar = ncOUT.createVariable('reO3h'        ,'d',('time','z','y','x')  ); ncvar[:] = D3

    ncvar = ncOUT.createVariable('index'        ,'i',('time','z','y','x')  ); ncvar[:] = index;
    ncvar = ncOUT.createVariable('index_inv'    ,'i',('waterpoints','dim3')); ncvar[:] = index_inv;

    ncOUT.close()

    os.system("mkdir -p " + test['Dir'] + '/BC/')
# Atmosphere NUT

    ATM_DATE=file2stringlist('KB/atm_date')

    for date in ATM_DATE:
        outfile = test['Dir'] + '/BC/ATM_' + date + '.nc'
        ncOUT   = NC.netcdf_file(outfile,'w')


        ncOUT.createDimension('lon'    ,jpi);
        ncOUT.createDimension('lat'    ,jpj);        
        ncvar = ncOUT.createVariable('atm_N1p'      ,'f',('lat','lon')                   ); ncvar[:] = 3.75866672509673e-09;
        ncvar = ncOUT.createVariable('atm_N3n'      ,'f',('lat','lon')                   ); ncvar[:] = 2.24183651189621e-07;

        ncOUT.close()

# Atmosphere CO2

    CO2_DATE=file2stringlist('KB/co2_date')


    D2=np.ones((jpj,jpi),float)


    for date in CO2_DATE:
        # Create CO2 file
        outfile = test['Dir'] + '/BC/CO2_' + date + '.nc'
        ncOUT   = NC.netcdf_file(outfile,'w')

        ncOUT.createDimension('lon',jpi);
        ncOUT.createDimension('lat',jpj);

        ncvar = ncOUT.createVariable('CO2'     ,'f',('lat','lon')                   ); ncvar[:] = 390.;

        ncOUT.close()

# TIN
    riv_N1p=np.zeros((jpj,jpi),dtype=float)-1.0; riv_N1p[0,0]=0.35*10**(-5)
    riv_N3n=np.zeros((jpj,jpi),dtype=float)-1.0; riv_N3n[0,0]=0.2*10**(-3)
    riv_N5s=np.zeros((jpj,jpi),dtype=float)-1.0; riv_N5s[0,0]=1.0*10**(-4)
    riv_O3c=np.zeros((jpj,jpi),dtype=float)-1.0; riv_O3c[0,0]=0.35
    riv_O3h=np.zeros((jpj,jpi),dtype=float)-1.0; riv_O3h[0,0]=0.01
    riv_O2o=np.zeros((jpj,jpi),dtype=float)-1.0; riv_O2o[0,0]=0.002

    for date in TIN_DATE:
        # Create RIV file
        outfile = test['Dir'] + '/BC/TIN_' + date + '.nc'
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
    for date in GIB_DATE:
        outfile = test['Dir'] + '/BC/GIB_' + date + '.nc'
        ncOUT   = NC.netcdf_file(outfile,'w')

        ncOUT.createDimension('lon'    ,jpi)
        ncOUT.createDimension('lat'    ,jpj)
        ncOUT.createDimension('dep'    ,jpk)

        N1p = np.ones((jpk,jpj,jpi),np.float32)*(-1.0)
        N3n = np.ones((jpk,jpj,jpi),np.float32)*(-1.0)
        N1p[:,:,0:first_quarter]=0.14
        N3n[:,:,0:first_quarter]=1.0
        ncvar = ncOUT.createVariable('gib_N1p','f',("dep","lat","lon") ); ncvar[:] = N1p
        ncvar = ncOUT.createVariable('gib_N3n','f',("dep","lat","lon") ); ncvar[:] = N3n
        ncOUT.close()
# Open Boundary
        BOUNDARY_CONCENTRATION={}
        BOUNDARY_CONCENTRATION['N5s'] = 2.0
        BOUNDARY_CONCENTRATION['N1p'] = 0.065 # mmol/m3
        BOUNDARY_CONCENTRATION['N3n']=   1.3 # mol/m3
        BOUNDARY_CONCENTRATION['O3c']= 28700 # mg/m3
        BOUNDARY_CONCENTRATION['O3h']=  2800 # mmol/m3
        OPEN_BOUNDARY=np.zeros((jpk,jpj,jpi),np.float32)
        I = jpi-2
        j_min=1
        j_max=int(jpj/2)
        for k in range(jpk):
            OPEN_BOUNDARY[k,j_min:j_max,I] = 0.1
        outfile = test['Dir'] + '/BC/OPE_yyyy0630-00:00:00.nc'
        ncOUT   = NC.netcdf_file(outfile,'w')
        ncOUT.createDimension('lon'    ,jpi);
        ncOUT.createDimension('lat'    ,jpj);
        ncOUT.createDimension('z'     ,jpk );
        for var in BOUNDARY_CONCENTRATION.keys():
            for k in range(jpk):
                OPEN_BOUNDARY[k,j_min:j_max,I] = BOUNDARY_CONCENTRATION[var]
            ncvar = ncOUT.createVariable(var,'f',('z','lat','lon')  )
            ncvar[:] = OPEN_BOUNDARY

        ncOUT.close()

