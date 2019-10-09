#prepare model and argo output for .nc files and matchup analysis
import netCDF4 as NC4
import numpy as np
 
def save_matchup(ncfile, PresCHL, Ed380_float, Ed412_float, Ed490_float, Ed380_model, Ed412_model, Ed490_model, timestr):

    modelfile = 'MATCHUP/' + ncfile
    ncmodel   = NC4.Dataset(modelfile,"w");
                
    ncmodel.createDimension('depth',     len(PresCHL));
    ncmodel.createDimension('wavelength', 3);
    
    floatstack = np.vstack((Ed380_float, Ed412_float, Ed490_float)).T
    modelstack = np.vstack((Ed380_model, Ed412_model, Ed490_model)).T
    
    ncvar = ncmodel.createVariable('Ed_float', 'f', ('depth', 'wavelength')); ncvar[:] = floatstack
    ncvar = ncmodel.createVariable('Ed_model', 'f', ('depth', 'wavelength')); ncvar[:] = modelstack
    
    setattr(ncmodel.variables,  'missing_value',-1.0 );     
    setattr(ncmodel.variables,  'long_name',  'Downward irradiance ' );     
    setattr(ncmodel.variables,  'unit',  '[uW/cm2/nm]' );  
    setattr(ncmodel.variables, 'time', timestr)   ;
        
    ncmodel.close()

    return ncmodel