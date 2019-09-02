import sys
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
    print(frame_start)

    frame_end   = frame_start + 4
    print(frame_end)

    
    # file to create cld data netcdf input file with original data converted to netcdf
    
    # read ECMWF ERA INTERIM
    
    indata_ECMWF_DIR = '/gpfs/scratch/userexternal/plazzari/BIOPTIMOD_HF/ECMWF_DATA/'
    
    filein_ecmwf= indata_ECMWF_DIR + 'ERAINTERIM_' + yyyymmdd[0:4] + '.nc'
    
    ncin=NC4.Dataset(filein_ecmwf,"r")
    
    # Mean sea level pressure
    # atmospheric pressure required by OASIM
    msl =  ncin.variables['msl'][frame_start:frame_end,1:181,:] # Mean sea level pressurea Pascal
    sp  =  ncin.variables['sp'][frame_start:frame_end,1:181,:] # surface pressure Pascal
    
    slporg=np.zeros(msl.shape)
    sporg=np.zeros(sp.shape)
    
    nTimes, nLat, nLon  = slporg.shape
    
    for t in range(nTimes):
        for j in range(nLat):
            for i  in range(nLon):
                 slporg[t,j,i] = msl[t,nLat-1-j,i] / 100. # convert to mb
                 sporg[t,j,i]  = sp[t,nLat-1-j,i]  / 100. # convert to mb
    
    
    set_nan_values(slporg, -1.0)
    set_nan_values(sporg, -1.0)
    
    # 10 metres wind speed
    u10 =  ncin.variables['u10'][frame_start:frame_end,1:181,:] # 10 metres u wind speed
    v10 =  ncin.variables['v10'][frame_start:frame_end,1:181,:] # 10 metres v wind speed

    si10 = np.sqrt(u10*u10 + v10*v10) 
    wsmorg=np.zeros(si10.shape)
    
    nTimes, nLat, nLon  = wsmorg.shape
    
    for t in range(nTimes):
        for j in range(nLat):
            for i  in range(nLon):
                 wsmorg[t,j,i] = si10[t,nLat-1-j,i] 
    
    set_nan_values(wsmorg, -1.0)
    
    # relative humidity and
    # precipitable water (absorption by water vapour)
    t2m  =  ncin.variables['t2m'][frame_start:frame_end,1:181,:] # 2m temperature K 
    d2m =  ncin.variables['d2m'][frame_start:frame_end,1:181,:] # dew temperature K
    
    a1 = 611.21 # Pascal
    a3 = 17.502 # dimensionless
    a4 = 32.19  # Kelvin
    To = 273.16 # Kelvin
    b1 = 0.14 * 0.01  # cm/Pascal
    b2 = 0.21   # cm
    
    rhorg=np.zeros(si10.shape)
    wvorg=np.zeros(si10.shape)
    
    nTimes, nLat, nLon  = rhorg.shape
    
    for t in range(nTimes):
        for j in range(nLat):
            for i  in range(nLon):
                 T            = t2m[t,nLat-1-j,i] 
                 Td           = d2m[t,nLat-1-j,i] 
                 es_Td        = a1 * np.exp( a3 * (Td-To)/(Td-a4) )
                 es_T         = a1 * np.exp( a3 *  (T-To)/(T-a4) )
                 rhorg[t,j,i] = 100. * es_Td /es_T
                 wvorg[t,j,i] = b1 * es_Td * sporg[t,j,i]/slporg[t,j,i] + b2
    
    set_nan_values(rhorg, -1.0)
    set_nan_values(wvorg, -1.0)
    #
    # computation of precipitable water [absorption by water vapour] from saturation water vapour pressure es_T as in Garrison and Adler (1990) and then used in Gregg and Casey (1990)
    #
    # note that saturation water pressure is computed in 3 ways:
    #
    # 1) in Garrison and Adler (1990) by the Tabata relation (Tabata, 1973)
    #    es_T = 10^(8.42926609 - 1827.17843/T - 71208.271/T^2) in mb
    #
    # 2) in Gregg and Carder (1990) they use the method of Lowe (Lowe, 1977)
    #    es_T = 1013.25*exp(13.3185*t - 1.9760 *t^2 - 0.6445*t^3 - 0.1299*t^4) in mb, where t=1-373.16/T
    #
    # 3) in ECMWF model (https://www.ecmwf.int/sites/default/files/elibrary/2015/9211-part-iv-physical-processes.pdf) by the Teten's formula (in Pa, 1 Pa = 0.01 mb) which is the one then used to compute rhorg above
    #
    # however, the three relations give the same results, as shown in saturation-water-vapour-pressure.png
    #
    # !!!! OCCORRE IMPORTARE ANCHE LA SURFACE PRESSURE sfcpr da ECMWF
    # !!!! CONTROLLARE SE SIA NECESSARIO DIVIDERE PER 100 IN QUANTO in GarAdl90 si definisce
    # e = h * es dove e=vapour pressure (ovvero es_Td) e h=relative humidity
    #
    
    # ozone
    
    ozo =  ncin.variables['tco3'][frame_start:frame_end,1:181,:] # total column ozone
    
    ozorg=np.zeros(si10.shape)
    nTimes, nLat, nLon  = ozorg.shape
    
    for t in range(nTimes):
        for j in range(nLat):
            for i  in range(nLon):
                 ozorg[t,j,i] = ozo[t,nLat-1-j,i] / 2.1414 / 0.00001
    
    set_nan_values(ozorg, -1.0)
    
    # write cloud file 
    
    outfile = 'output/opt' + yyyymmdd  + '_ECMWF.nc'
    ncOUT   = NC4.Dataset(outfile,"w");
            
    ncOUT.createDimension('lon',nLon);
    ncOUT.createDimension('lat',nLat);
    ncOUT.createDimension('time',nTimes);
    
#   ncvar = ncOUT.createVariable('slporg','f',('time','lat','lon')); ncvar[:] = sporg; 
    sporg_ave = np.mean(sporg,axis=0)
    ncvar = ncOUT.createVariable('slporg','f',('lat','lon')); ncvar[:] = sporg_ave; 
    setattr(ncOUT.variables['slporg'], 'missing_value',-1.0 );
    setattr(ncOUT.variables['slporg'], 'input_file', filein_ecmwf );
    setattr(ncOUT.variables['slporg'], 'long_name',  'surface pressure' );
    setattr(ncOUT.variables['slporg'], 'unit',  '[mb]' );
    
#   ncvar = ncOUT.createVariable('wsmorg','f',('time','lat','lon')); ncvar[:] = wsmorg; 
    wsmorg_ave = np.mean(wsmorg,axis=0)
    ncvar = ncOUT.createVariable('wsmorg','f',('lat','lon')); ncvar[:] = wsmorg_ave; 
    setattr(ncOUT.variables['wsmorg'], 'missing_value',-1.0 );     
    setattr(ncOUT.variables['wsmorg'], 'input_file', filein_ecmwf );
    setattr(ncOUT.variables['wsmorg'], 'long_name',  '10 m wind speed' );
    setattr(ncOUT.variables['wsmorg'], 'unit',  '[m.s-1]' );
    
#   ncvar = ncOUT.createVariable('rhorg','f',('time','lat','lon')); ncvar[:] = rhorg; 
    rhorg_ave = np.mean(rhorg,axis=0)
    ncvar = ncOUT.createVariable('rhorg','f',('lat','lon')); ncvar[:] = rhorg_ave; 
    setattr(ncOUT.variables['rhorg'], 'missing_value',-1.0 );     
    setattr(ncOUT.variables['rhorg'], 'input_file', filein_ecmwf );
    setattr(ncOUT.variables['rhorg'], 'long_name',  'relative humidity');
    setattr(ncOUT.variables['rhorg'], 'unit',  '[%]' );
    
#   ncvar = ncOUT.createVariable('ozorg' ,'f',('time','lat','lon')); ncvar[:] = ozorg; 
    ozorg_ave = np.mean(ozorg,axis=0)
    ncvar = ncOUT.createVariable('ozorg' ,'f',('lat','lon')); ncvar[:] = ozorg_ave; 
    setattr(ncOUT.variables['ozorg'], 'missing_value',-1.0 );     
    setattr(ncOUT.variables['ozorg'], 'input_file', filein_ecmwf );
    setattr(ncOUT.variables['ozorg'], 'long_name',  'ozone');
    setattr(ncOUT.variables['ozorg'], 'unit',  '[DU]' );
    
#   ncvar = ncOUT.createVariable('wvorg','f',('time','lat','lon')); ncvar[:] = wvorg; 
    wvorg_ave = np.mean(wvorg,axis=0)
    ncvar = ncOUT.createVariable('wvorg','f',('lat','lon')); ncvar[:] = wvorg_ave; 
    setattr(ncOUT.variables['wvorg'], 'missing_value',-1.0 );     
    setattr(ncOUT.variables['wvorg'], 'input_file', filein_ecmwf );
    setattr(ncOUT.variables['wvorg'], 'long_name',  'precipitable water');
    setattr(ncOUT.variables['wvorg'], 'unit',  '[cm]' );
    
    ncOUT.close()

if __name__ == "__main__":
    main()

