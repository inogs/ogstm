from basins import V2 as OGS
from commons.layer import Layer
from commons.time_interval import TimeInterval
from instruments import optbio_float_2019
from instruments import var_conversions
from matchup_manager import Matchup_Manager
import numpy as np
import os
from profiler import *
import scipy.io.netcdf as NC
import netCDF4 as NC4

def PFT_calc(CHL):
    PFT_1 = 0.25*CHL
    PFT_2 = 0.25*CHL
    PFT_3 = 0.25*CHL
    PFT_4 = 0.25*CHL
    return PFT_1, PFT_2, PFT_3, PFT_4

def NAP_calc(CHL):
    NAP = 0.*CHL
    return NAP

def CDOM_calc(CHL):
    CDOM = 0.*CHL
    return CDOM

wl =np.array( [250., 325., 350., 375., 400.,   425.,  450.,  475.,  500.,  525.,  550.,  575.,  600.,  625.,  650.,  675., 700.,
               725., 775., 850., 950., 1050., 1150., 1250., 1350., 1450., 1550., 1650., 1750., 1900., 2200., 2900., 3700.])
wl_int = wl.astype(np.int64)

#Indices for float wavelengths (380, 412, 490) will be 3, 4 and 7 respectively

wl_str = [np.str(wl_int[i]) for i in range(len(wl_int))] 
str_Ed = ['Ed_' + np.str(wl_str[i]) for i in range(len(wl_str))]
str_Es = ['Es_' + np.str(wl_str[i]) for i in range(len(wl_str))]

maskfile    = '/galileo/home/userexternal/eterzic0/OASIM_POSTPROC/ARGO_MATCHUP/CODES/meshmask.nc'
ncIN=NC.netcdf_file(maskfile,"r")
nav_lev = ncIN.variables['nav_lev'].data.copy()
ncIN.close()

M = Matchup_Manager(ALL_PROFILES,TL,BASEDIR)

TI = TimeInterval("20120101", "20171231","%Y%m%d")
variable='P_l'
varname=var_conversions.FLOAT_OPT_VARS_2019[variable]

Profilelist=optbio_float_2019.FloatSelector(varname,TI , OGS.med)
p = Profilelist[0]

#Read OASIM Ed and Es data
List_Ed = [M.getMatchups([p], nav_lev, modelvar).subset(Layer(0,1.5)) for modelvar in str_Ed]
List_Es = [M.getMatchups([p], nav_lev, modelvar).subset(Layer(0,1.5)) for modelvar in str_Es]

Ed = [List_Ed[i].Model[0] for i in range(len(List_Ed))]
Es = [List_Es[i].Model[0] for i in range(len(List_Es))]

np.savetxt('OASIM.txt', np.c_[Ed, Es])

#Read BGC-ARGO profiles
PresCHL, CHLz,    Qc = p.read('CHL')
Pres380, Ed_380,  Qc = p.read('IRR_380')
Pres412, Ed_412,  Qc = p.read('IRR_412')
Pres490, Ed_490,  Qc = p.read('IRR_490')
PresPAR, PAR,     Qc = p.read('PAR')

Lon = p.lon
Lat = p.lat

timestr = p.time.strftime("%Y%m%d-%H:%M:%S")

nLevels = len(PresCHL)

init_rows = str(timestr) + '\n' + str(Lat) + '\n' + str(nLevels)

# Calculate PFT
PFT1, PFT2, PFT3, PFT4 = PFT_calc(CHLz)

# Calculate NAP and CDOM
NAP  = NAP_calc(CHLz)
CDOM = CDOM_calc(CHLz)

file_cols = np.vstack((PresCHL, PFT1, PFT2, PFT3, PFT4, CDOM, NAP)).T
np.savetxt('IOP.txt', file_cols, header = init_rows, delimiter='\t', comments='')

if len(PresCHL) < 15:
    badstr = '_BAD'
    floatname = p.ID() + '_BAD.nc'
else:
    floatname = p.ID() + '.nc'

command='./compute.xx ' + 'OASIM.txt ' + 'IOP.txt ' + str(floatname) 
os.system(command)

# Move them to a new directory

#movefiles = 'mv *.nc NCOUT/'
#os.system(movefiles)

ncin=NC4.Dataset(floatname,"r")

Ed380 =  np.array(ncin.variables['Edz'][3,1:] + ncin.variables['Esz'][3,1:])  
Ed412 =  np.array(ncin.variables['Edz'][4,1:] + ncin.variables['Esz'][4,1:]) 
Ed490 =  np.array(ncin.variables['Edz'][7,1:] + ncin.variables['Esz'][7,1:]) 

Ed380_model = Ed380 * 4  # = 10**(-6) / (10**(-4) * 25) 
Ed412_model = Ed412 * 4  #  W/m2 to muW/cm2
Ed490_model = Ed490 * 4  

#Interpolate Ed380 on CHL (OASIM model) depth quotes

Ed380_float = np.interp(PresCHL, Pres380, Ed_380)
Ed412_float = np.interp(PresCHL, Pres412, Ed_412)
Ed490_float = np.interp(PresCHL, Pres490, Ed_490)

# Save irradiance profiles of floats and the ones of models at the corresponding wavelengths

#Save Es
outfile = 'MATCHUP/' + floatname 
ncOUT   = NC4.Dataset(outfile,"w");
            
ncOUT.createDimension('depth',     len(PresCHL));
ncOUT.createDimension('wavelength',3);

ncvar = ncOUT.createVariable('Ed_float', 'f', ('depth', 'wavelength')); ncvar[:] = [Ed380_float, Ed412_float, Ed490_float]
ncvar = ncOUT.createVariable('Ed_model', 'f', ('depth', 'wavelength')); ncvar[:] = [Ed380_model, Ed412_model, Ed490_model]

setattr(ncOUT.variables,  'missing_value',-1.0 );     
setattr(ncOUT.variables,  'long_name',  'Downward irradiance ' );     
setattr(ncOUT.variables,  'unit',  '[uW/cm2/nm]' );     
    
ncOUT.close()
    
# Copy the model output to a separate directory

movefiles = 'cp ' + str(floatname) + ' NCOUT/'
os.system(movefiles)


