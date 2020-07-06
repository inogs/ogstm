#!/bin/env python

import netCDF4 as NC4
import numpy as np

''' Here you put all the functions, also for the IOPs , T-S corrections, etc. '''

def findVars(Varlist):
    allvars=[' CHL', 'IRR_380', 'IRR_412', 'IRR_490', 'PAR']
    if len(Varlist)==0: return False    
    
    for var in allvars:
    	if not var in Varlist:
    		return False
    return True

'''
Prepare model and BGC-Argo output for .nc files and match-up analysis
'''

def save_matchup(ncfile, PresCHL, Ed380_float, Ed412_float, Ed490_float, Ed380_model, Ed412_model, Ed490_model, timestr):

    modelfile = 'MATCHUP/' + ncfile
    ncmodel   = NC4.Dataset(modelfile,"w");
                
    ncdepth = ncmodel.createDimension('depth',     len(PresCHL));
    ncwave  = ncmodel.createDimension('wavelength', 3);

    setattr(ncmodel, 'time', timestr);
    
    ncDepth = ncmodel.createVariable('depth', 'f', ('depth')); 
    setattr(ncDepth, 'unit',  '[m]' );
    ncDepth[:] = PresCHL
    
    ncEdf = ncmodel.createVariable('Ed_float', 'f', ('depth', 'wavelength'));
    setattr(ncEdf, 'missing_value',-1.0 );     
    setattr(ncEdf, 'long_name',  'Downward irradiance ' );     
    setattr(ncEdf, 'unit',  '[uW/cm2/nm]' );

    ncEdf[:] = np.vstack((Ed380_float, Ed412_float, Ed490_float)).T
    
    ncEdm = ncmodel.createVariable('Ed_model', 'f', ('depth', 'wavelength'));
    setattr(ncEdm, 'missing_value',-1.0 );     
    setattr(ncEdm, 'long_name',  'Downward irradiance ' );     
    setattr(ncEdm, 'unit',  '[uW/cm2/nm]' );

    ncEdm[:] = np.vstack((Ed380_model, Ed412_model, Ed490_model)).T
    
    ncmodel.close()

    return ncmodel