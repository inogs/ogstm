import numpy as np
from mydtype import *
import scipy.io.netcdf as NC
import os


def create_grid(test):
    OUTDIR=test['Dir'] + "/DA_static_data/3DVAR/GRID/"
    filename = OUTDIR +  "grid1.nc"
    jpi=test['jpi'];
    jpj=test['jpj'];
    jpk=test['jpk'];
    
    maskfile=test['Dir'] + '/meshmask.nc'

    M=NC.netcdf_file(maskfile,"r")

    nav_lon     =  M.variables['nav_lon'].data.copy()
    nav_lat     =  M.variables['nav_lat'].data[0,0,:,:].copy()
    gdept   =  M.variables['gdept'].data[0,:,0,0].copy()
    gdepw   =  M.variables['gdepw'].data[0,:,0,0].copy()
    e1v     =  M.variables['e1v'  ].data[0,0,:,:].copy()
    e2u     =  M.variables['e2u'  ].data[0,0,:,:].copy()
    e3t     =  M.variables['e3t'  ].data[0,:,0,0].copy()
    
    M.close()
    
    
    
    NCout = NC.netcdf_file(filename,"w")
    NCout.createDimension("im",jpi)
    NCout.createDimension("jm",jpj)
    NCout.createDimension("km",jpk)
    
    ncvar = NCout.createVariable("lon", 'f', ('im', 'jm'))
     
    NCout.createVariable("lon", 'f', ('im', 'jm'))
    NCout.createVariable("lon", 'f', ('im', 'jm'))
    
    
    NCout.close()

    


def create_variance_for_misfit(test):
    jpi=test['jpi']
    jpj=test['jpj']
    dx=test['dx']
    dy=test['dy']
    
    lon = np.zeros((jpi,), np.float32)
    lat = np.zeros((jpj,), np.float32)
    
    for ji in range(jpi): lon[ji]=test['lon0']-jpi/2.*dx+dx/2. + dx*ji;
    for jj in range(jpj): lat[jj]=test['lat0']-jpi/2.*dy+dy/2. + dy*jj;
   
    OUTDIR=test['Dir'] + "/DA_static_data/MISFIT/VAR2D/"
    os.system("mkdir -p " + OUTDIR)
    
    for month in range(12):
        filename = OUTDIR +  "/var2D.%02d.nc"  %(month+1)
    
        NCout = NC.netcdf_file(filename,"w")
        NCout.createDimension("lon", jpi)
        NCout.createDimension("lat", jpj)
        
        ncvar=NCout.createVariable("variance", 'f', ('lat','lon'))
        setattr(ncvar,"missing_value", 1.e+20)
        ncvar[:]= 0.0002 * (np.random.rand(jpj,jpi).astype(np.float32) -0.5) + 0.0005
                
        NCout.close()


def create_eofs_for_3dvar(test):
    nreg = 14
    neof = 4
    nlev = 26 #test['jpk']
    
    EVA = 0.08 * (np.random.rand(neof,nreg).astype(np.float32) - 0.5 ) + 0.15
    EVC = 12 * np.random.rand(neof,nlev,nreg).astype(np.float32) -6
        
    OUTDIR=test['Dir'] + "/DA_static_data/3D_VAR/EOF/"
    os.system("mkdir -p " + OUTDIR)
    for month in range(12):
        filename = OUTDIR +  "eof.%02d.nc"  %(month+1)
        NCout = NC.netcdf_file(filename,"w")
        NCout.createDimension("nreg", nreg)
        NCout.createDimension("nlev", nlev)
        NCout.createDimension("neof", neof)
        
        ncvar = NCout.createVariable("eva", 'f', ('neof','nreg'))
        ncvar[:] = EVA
        ncvar = NCout.createVariable("evc", 'f', ('neof','nlev','nreg'))
        ncvar[:] = EVC                
        NCout.close()
    OUTDIR=test['Dir'] + "/DA_static_data/3D_VAR/GRID/"
    os.system("mkdir -p " + OUTDIR)    

def create_misfit_for_3dvar(test): 
    jpi=test['jpi']
    jpj=test['jpj']
    dx=test['dx']
    dy=test['dy']
    
    lon = np.zeros((jpi,), np.float32)
    lat = np.zeros((jpj,), np.float32)
    
    for ji in range(jpi): lon[ji]=test['lon0']-jpi/2.*dx+dx/2. + dx*ji;
    for jj in range(jpj): lat[jj]=test['lat0']-jpi/2.*dy+dy/2. + dy*jj;
    
    MIS = 0.20 * (np.random.rand(jpj,jpi).astype(np.float32) -0.5) + 0.25
    ERR = 0.20 * (np.random.rand(jpj,jpi).astype(np.float32) -0.5) + 0.45
    NCout = NC.netcdf_file(test['Dir']+ '/chl_misfit.nc',"w")
    NCout.createDimension("lon", jpi)
    NCout.createDimension("lat", jpj)
    NCout.createDimension("time",  1)
    NCout.createDimension("depth", 1)
    
    ncvar=NCout.createVariable("time", 'f', ('time',))
    ncvar[:] = 1
    ncvar=NCout.createVariable("depth", 'f', ('depth',))
    ncvar[:] = 0

    
    ncvar=NCout.createVariable("lat", 'f', ('lat',))
    ncvar[:] = lat
    ncvar=NCout.createVariable("lon", 'f', ('lon',))
    ncvar[:] = lon
    
    ncvar=NCout.createVariable("misfchl", 'f', ('lat','lon'))
    setattr(ncvar,"missing_value", 1.e+20)
    ncvar[:]= MIS
    ncvar=NCout.createVariable("errchl", 'f', ('lat','lon'))
    setattr(ncvar,"missing_value", 1.e+20)
    ncvar[:]= ERR    
    
    
    NCout.close()          
    
    
    NCout.close()
    
def create_sat_files(test):
    OUTDIR=test['Dir'] + "/SATELLITE/"
    os.system("mkdir -p " + OUTDIR)

    jpi=test['jpi']
    jpj=test['jpj']
    dx=test['dx']
    dy=test['dy']
    
    lon = np.zeros((jpi,), np.float32)
    lat = np.zeros((jpj,), np.float32)
    
    for ji in range(jpi): lon[ji]=test['lon0']-jpi/2.*dx+dx/2. + dx*ji;
    for jj in range(jpj): lat[jj]=test['lat0']-jpi/2.*dy+dy/2. + dy*jj;
    
    SAT = 0.20 * (np.random.rand(jpj,jpi).astype(np.float32) -0.5) + 0.25

    DATElist = ["20000102", "20000105", "20000109"]
    
    for date in DATElist:
        filename = OUTDIR +  date + "_d-OC_CNR-L4-CHL-MedOC4_SAM_7KM-MED-REP-v02.nc"
    
        NCout = NC.netcdf_file(filename,"w")
        NCout.createDimension("lon", jpi)
        NCout.createDimension("lat", jpj)
        NCout.createDimension("time",  1)
        
        ncvar=NCout.createVariable("lat", 'f', ('lat',))
        ncvar[:] = lat
        ncvar=NCout.createVariable("lon", 'f', ('lon',))
        ncvar[:] = lon
        
        ncvar=NCout.createVariable("lchlm", 'f', ('lat','lon'))
        setattr(ncvar,"missing_value", 1.e+20)
        ncvar[:]= SAT
        NCout.close()       
    

def create_dataset(test):
    create_sat_files(test)
    create_variance_for_misfit(test)
    create_eofs_for_3dvar(test)
    create_misfit_for_3dvar(test)
    print "For the mesh lauch this command:"
    print "./createGridDA meshmask.nc submask.nc ESO -5.7 BFM_grid.nc"
    print "ncks -d km,1,26 BFM_grid.nc -O DA_static_data/3D_VAR/GRID/BFM_grid.nc"
    print "make sure to have nreg=14 and neof =4 in var_3d_nml"
    
        
