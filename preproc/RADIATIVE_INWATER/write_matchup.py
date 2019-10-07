#prepare model and argo output for .nc files and matchup analysis

def write_matchup(ncfile):
    #Indices for float wavelengths (380, 412, 490) will be 3, 4 and 7 respectively
    
    ncin=NC4.Dataset(ncfile,"r")
    
    Ed380_model  =  np.array(ncin.variables['Edz'][3,1:] + ncin.variables['Esz'][3,1:])  * 4 # = 10**(-6) / (10**(-4) * 25) 
    Ed412_model  =  np.array(ncin.variables['Edz'][4,1:] + ncin.variables['Esz'][4,1:])  * 4 #  W/m2 to muW/cm2
    Ed490_model  =  np.array(ncin.variables['Edz'][7,1:] + ncin.variables['Esz'][7,1:])  * 4
    
    ncin.close()
    #Interpolate Ed380 on CHL (OASIM model) depth quotes
    
    Ed380_float = np.interp(PresCHL, Pres380, Ed_380)
    Ed412_float = np.interp(PresCHL, Pres412, Ed_412)
    Ed490_float = np.interp(PresCHL, Pres490, Ed_490)

    return Ed380_float, Ed380_model, Ed412_float, Ed412_model, Ed490_float, Ed490_model
 
def save_matchup(ncfile):

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
        
    ncmodel.close()

    return ncmodel