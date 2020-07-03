import argparse
import datetime
import netCDF4 as NC4
import numpy as np
import os.path
import sys

def argument():
    parser = argparse.ArgumentParser(description = '''
    Creates density plots of surface irradiance data 
    from BGC-Argo floats and the atmospheric OASIM output
    ''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(   '--date', '-d',
                                type = str,
                                required = True,
                                help = ''' Specify the date for which to perform the aveScan (yyyymmdd format).'''
                                )

    return parser.parse_args()

args = argument()

wl =np.array( [250., 325., 350., 375., 400.,   425.,  450.,  475.,  500.,  525.,  550.,  575.,  600.,  625.,  650.,  675., 700.,
               725., 775., 850., 950., 1050., 1150., 1250., 1350., 1450., 1550., 1650., 1750., 1900., 2200., 2900., 3700.])

wl_int = wl.astype(np.int64)

wl_str = [np.str(wl_int[i]) for i in range(len(wl_int))] 
    
str_Ed = ['Ed_' + np.str(wl_str[i]) for i in range(len(wl_str))]
str_Es = ['Es_' + np.str(wl_str[i]) for i in range(len(wl_str))]

filein = '/gpfs/scratch/userexternal/eterzic0/OASIM_HF_INWATER/DATA/' + 'rad_0m' + args.date + '.nc' 
if os.path.isfile(filein):
    print([ 'reading model file ' + filein ])
else: 
    print([ 'reading model file ' + filein + ' --> not found' ])

ncin=NC4.Dataset(filein,"r")


for iwl in range(len(wl)):   # Loop over OASIM wavelengths

    Ed =  (ncin.variables['Ed_0m'][:,iwl,:,:]) 
    Es =  (ncin.variables['Es_0m'][:,iwl,:,:]) 
    
    if os.path.isfile(filein):
        print([ 'reading model file ' + filein ])
    else:
        print([ 'reading model file ' + filein + ' --> not found' ])
    
    jpi=360
    jpj=180 
    
    # Read year, month, date from the args.date object
    date_object = datetime.datetime.strptime(args.date, "%Y%m%d")
    step = datetime.timedelta(minutes=15)
    d1 = date_object

    for itime in range(96):  # Loop over time - 15-min intervals

        d1 = d1 + step
        HHMMSS = d1.strftime('%H:%M:%S')

        if itime > 38 and itime < 56:   # indices between 10 and 14:00 local time
            print(itime, HHMMSS) 
        else:
            continue

        # Save Ed
        outfile = '/gpfs/scratch/userexternal/eterzic0/OASIM_HF_INWATER/AVEDATA/ave.' + args.date  +  '-' + HHMMSS + '.Ed_' + str(wl_int[iwl]) + '.nc'
        ncOUT   = NC4.Dataset(outfile,"w");

        ncOUT.createDimension('lon',jpi);
        ncOUT.createDimension('lat',jpj);
        ncOUT.createDimension('time',1);
        
        ncvar = ncOUT.createVariable(str_Ed[iwl],'f',('time','lat','lon')); ncvar[:] = Ed[itime,:,:]
        setattr(ncOUT.variables[str_Ed[iwl]],  'missing_value',-1.0 );     
        setattr(ncOUT.variables[str_Ed[iwl]],  'long_name',  'Downward irradiance ' );     
        setattr(ncOUT.variables[str_Ed[iwl]],  'unit',  '[uW/cm2/nm]' );     
        
        ncOUT.close()
        
        # Save Es
        outfile = '/gpfs/scratch/userexternal/eterzic0/OASIM_HF_INWATER/AVEDATA/ave.' + args.date   +  '-' + HHMMSS + '.Es_' + str(wl_int[iwl]) + '.nc'
        ncOUT   = NC4.Dataset(outfile,"w");
                
        ncOUT.createDimension('lon',jpi);
        ncOUT.createDimension('lat',jpj);
        ncOUT.createDimension('time',1);
        
        ncvar = ncOUT.createVariable(str_Es[iwl],'f',('time','lat','lon')); ncvar[:] = Es[itime,:,:]
        setattr(ncOUT.variables[str_Es[iwl]],  'missing_value',-1.0 );     
        setattr(ncOUT.variables[str_Es[iwl]],  'long_name',  'Downward irradiance ' );     
        setattr(ncOUT.variables[str_Es[iwl]],  'unit',  '[uW/cm2/nm]' );     
        
        ncOUT.close()
        
