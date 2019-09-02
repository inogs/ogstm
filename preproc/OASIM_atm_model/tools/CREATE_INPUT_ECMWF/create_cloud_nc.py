import sys
import os.path
from datetime import datetime
import netCDF4 as NC4
import numpy as np
    
def set_nan_values(M, fill_val):
    
    nTimes, nLat, nLon  = M.shape
    
    for t in range(nTimes):
        for j in range(nLat):
            for i  in range(nLon):
                 if np.isnan(M[t,j,i]): 
                    M[t,j,i] = fill_val 
                 if M[t,j,i] < 0.: 
                    M[t,j,i] = fill_val 
    
def main():
    # print command line arguments
    for i,arg in enumerate(sys.argv[1:]):
        if i == 0:
           yyyymmdd = arg
        if i >  0:
           raise TypeError(" max #arg = 0") 


    yyyy = np.int(yyyymmdd[0:4])

    mm   = np.int(yyyymmdd[4:6]) 

    dd   = np.int(yyyymmdd[6:8]) 

    frame_start = ( datetime(yyyy, mm, dd).timetuple().tm_yday - 1 ) * 4

    frame_end   = frame_start + 4

    # file to create cld data netcdf input file with original data converted to netcdf
    
    # read MODIS data variables cldtcm and cdrem
    
    indata_MODIS_DIR = '/gpfs/scratch/userexternal/plazzari/BIOPTIMOD_HF/MODIS_DATA/'
    
    filein_modis= indata_MODIS_DIR + 'modcld' + yyyymmdd[0:4] + '.nc'

    if os.path.isfile(filein_modis): 
       print([ 'reading modis file ' + filein_modis ])
    else: 
       print([ 'reading modis file ' + filein_modis + ' --> not found' ])
       filein_modis = indata_MODIS_DIR + 'modcld0000.nc'
       print([ 'reading modis file ' + filein_modis ])

    ncin=NC4.Dataset(filein_modis,"r")
    cldtcm =  ncin.variables['cldtcm'][mm-1,:,:]
    cldtc=np.zeros((4,cldtcm.shape[0],cldtcm.shape[1]))
    nTimes, nLat, nLon  = cldtc.shape
    print(nTimes)
    print(nLat)
    print(nLon)
    
    for t in range(nTimes):
        for j in range(nLat):
            for i  in range(nLon):
                 cldtc[t,j,i] = cldtcm[j,i]
    
    set_nan_values(cldtc, -1.0)
                    
    
    cdrem = ncin.variables['cdrem'][mm-1,:,:]
    
    cdre = np.zeros((4,cdrem.shape[0],cdrem.shape[1]))

    nTimes, nLat, nLon  = cdre.shape
    print(nTimes)
    print(nLat)
    print(nLon)
    
    for t in range(nTimes):
        for j in range(nLat):
            for i  in range(nLon):
                 cdre[t,j,i] = cdrem[j,i]
    
    set_nan_values(cdre, -1.0)
    
    ncin.close()
    
    # read ECMWF ERA INTERIM
    
    indata_ECMWF_DIR = '/gpfs/scratch/userexternal/plazzari/BIOPTIMOD_HF/ECMWF_DATA/'
    
    filein_ecmwf= indata_ECMWF_DIR + 'ERAINTERIM_' + yyyymmdd[0:4] + '.nc'
    
    ncin=NC4.Dataset(filein_ecmwf,"r")
    
    ccove =  ncin.variables['tcc'][frame_start:frame_end,1:181,:]
    
    ccov=np.zeros(ccove.shape)
    
    nTimes, nLat, nLon  = ccove.shape
    
    for t in range(nTimes):
        for j in range(nLat):
            for i  in range(nLon):
                 ccov[t,j,i] = ccove[t,nLat-1-j,i] * 100.
    
    set_nan_values(ccov, -1.0)
    
    filein_ecmwf= indata_ECMWF_DIR + 'ERAINTERIM_' + yyyymmdd[0:4] + '.nc'
    
    ncin=NC4.Dataset(filein_ecmwf,"r")
    
    rlwpi =  ncin.variables['p56.162'][frame_start:frame_end,1:181,:] # vilw
#   rlwpi =  ncin.variables['tclw'][:,1:181,:] # not found in daily files
    
    rlwp=np.zeros(rlwpi.shape)
    nTimes, nLat, nLon  = rlwpi.shape
    
    for t in range(nTimes):
        for j in range(nLat):
            for i  in range(nLon):
                 rlwp[t,j,i] = rlwpi[t,j,i] * 1000.
    
    set_nan_values(rlwp, -99.0)
    
    ncin.close()
    
    # write cloud file 
    
    outfile = 'output/clouds' + yyyymmdd  + '_ECMWF.nc'
    ncOUT   = NC4.Dataset(outfile,"w");
            
    ncOUT.createDimension('lon',nLon);
    ncOUT.createDimension('lat',nLat);
    ncOUT.createDimension('time',nTimes);
    
#    ncvar = ncOUT.createVariable('ccov' ,'f',('time','lat','lon')); ncvar[:] = ccov; 
    ccov_ave=np.mean(ccov,axis=0)
    ncvar = ncOUT.createVariable('ccov' ,'f',('lat','lon')); ncvar[:] = ccov_ave; 
    setattr(ncOUT.variables['ccov'], 'missing_value',-1.0 );     
    setattr(ncOUT.variables['ccov'], 'input_file', filein_ecmwf );     
    setattr(ncOUT.variables['ccov'], 'long_name',  'cloud cover' );     
    setattr(ncOUT.variables['ccov'], 'unit',  '[%]' );     
    
#   ncvar = ncOUT.createVariable('cldtc','f',('time','lat','lon')); ncvar[:] = cldtc; 
    cldtc_ave = np.mean(cldtc,axis=0)
    ncvar = ncOUT.createVariable('cldtc','f',('lat','lon')); ncvar[:] = cldtc_ave; 
    setattr(ncOUT.variables['cldtc'], 'missing_value',-1.0 );     
    setattr(ncOUT.variables['cldtc'], 'input_file',  filein_modis );     
    setattr(ncOUT.variables['cldtc'], 'long_name',  'cloud optical thickness' );     
    setattr(ncOUT.variables['cldtc'], 'unit',  '[]' );     
    
#   ncvar = ncOUT.createVariable('rlwp','f',('time','lat','lon')); ncvar[:] = rlwp; 
    rlwp_ave=np.mean(rlwp,axis=0)
    ncvar = ncOUT.createVariable('rlwp','f',('lat','lon')); ncvar[:] = rlwp_ave; 
    setattr(ncOUT.variables['rlwp'], 'missing_value',-99.0 );     
    setattr(ncOUT.variables['rlwp'], 'input_file',  filein_ecmwf );     
    setattr(ncOUT.variables['rlwp'], 'long_name',  'cloud liquid water path' );     
    setattr(ncOUT.variables['rlwp'], 'unit',  '[g.m-2]' );     
    
    
#   ncvar = ncOUT.createVariable('cdre','f',('time','lat','lon')); ncvar[:] = cdre; 
    cdre_ave =np.mean(cdre,axis=0)
    ncvar = ncOUT.createVariable('cdre','f',('lat','lon')); ncvar[:] = cdre_ave; 
    setattr(ncOUT.variables['cdre'], 'missing_value',-1.0 );     
    setattr(ncOUT.variables['cdre'], 'input_file',  filein_modis );     
    setattr(ncOUT.variables['cdre'], 'long_name',  'cloud droplet effective radius' );     
    setattr(ncOUT.variables['cdre'], 'unit',  '[um]' );     
    
    ncOUT.close()

if __name__ == "__main__":
    main()

