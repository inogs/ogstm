import argparse
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

filein = '/gpfs/scratch/userexternal/eterzic0/OASIM/GMT_ALL/UNZIPPED/' + 'rad_0m' + args.date + '.nc' 
if os.path.isfile(filein):
    print([ 'reading model file ' + filein ])
else: 
    print([ 'reading model file ' + filein + ' --> not found' ])

ncin=NC4.Dataset(filein,"r")

for i in range(len(wl)):
    Ed =  (ncin.variables['Ed_0m'][5,i,:,:] + ncin.variables['Ed_0m'][6,i,:,:]) /2.0 # average over 10-14 day hours
    Es =  (ncin.variables['Es_0m'][5,i,:,:] + ncin.variables['Es_0m'][6,i,:,:]) /2.0 # average over 10-14 day hours
    
    if os.path.isfile(filein):
        print([ 'reading model file ' + filein ])
    else:
        print([ 'reading model file ' + filein + ' --> not found' ])
    
    jpi=360
    jpj=180 
    
    #Save Ed
    outfile = '/gpfs/scratch/userexternal/eterzic0/RADIATIVE_INWATER/INDATA/ave.' + args.date  + '-12:00:00.Ed_' + str(wl_int[i]) + '.nc'
    ncOUT   = NC4.Dataset(outfile,"w");
            
    ncOUT.createDimension('lon',jpi);
    ncOUT.createDimension('lat',jpj);
    ncOUT.createDimension('time',1);
    
    ncvar = ncOUT.createVariable(str_Ed[i],'f',('time','lat','lon')); ncvar[:] = Ed
    setattr(ncOUT.variables[str_Ed[i]],  'missing_value',-1.0 );     
    setattr(ncOUT.variables[str_Ed[i]],  'long_name',  'Downward irradiance ' );     
    setattr(ncOUT.variables[str_Ed[i]],  'unit',  '[uW/cm2/nm]' );     
    
    ncOUT.close()
    
    #Save Es
    outfile = '/gpfs/scratch/userexternal/eterzic0/RADIATIVE_INWATER/INDATA/ave.' + args.date  + '-12:00:00.Es_' + str(wl_int[i]) + '.nc'
    ncOUT   = NC4.Dataset(outfile,"w");
            
    ncOUT.createDimension('lon',jpi);
    ncOUT.createDimension('lat',jpj);
    ncOUT.createDimension('time',1);
    
    ncvar = ncOUT.createVariable(str_Es[i],'f',('time','lat','lon')); ncvar[:] = Es
    setattr(ncOUT.variables[str_Es[i]],  'missing_value',-1.0 );     
    setattr(ncOUT.variables[str_Es[i]],  'long_name',  'Downward irradiance ' );     
    setattr(ncOUT.variables[str_Es[i]],  'unit',  '[uW/cm2/nm]' );     
    
    ncOUT.close()
    
