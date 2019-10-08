from basins import V2 as OGS
from create_IOP import *
from commons.layer import Layer
from instruments import optbio_float_2019
from instruments import var_conversions
from instruments.matchup_manager import Matchup_Manager
from matchup import statistics
from model_params import *
import numpy as np
import os
from profiler import *
import scipy.io.netcdf as NC
from write_matchup import *


#######  Read model data that corresponds to the position and time of BGC-Argo floats ########

M = Matchup_Manager(ALL_PROFILES,TL,BASEDIR)

TI = TimeInterval("20120101", "20171231","%Y%m%d")
variable='P_l'
varname=var_conversions.FLOAT_OPT_VARS_2019[variable]

Profilelist=optbio_float_2019.FloatSelector(varname,TI , OGS.med)
p = Profilelist[15]
profile_ID = p.ID()


######################## phase 1. write OASIM.txt file #######################################

List_Ed = [M.getMatchups([p], nav_lev, modelvar).subset(Layer(0,1.5)) for modelvar in str_Ed]
List_Es = [M.getMatchups([p], nav_lev, modelvar).subset(Layer(0,1.5)) for modelvar in str_Es]

Ed = [List_Ed[i].Model[0] for i in range(len(List_Ed))]
Es = [List_Es[i].Model[0] for i in range(len(List_Es))]

np.savetxt(profile_ID + '_OASIM.txt', np.c_[Ed, Es])


###############################################################################################
####################### phase 2. Read BGC-ARGO profiles #######################################
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


################################################################################################
####################### phase 3. Calculate and save IOPs  ######################################
PFT1, PFT2, PFT3, PFT4 = PFT_calc(CHLz, 0.25, 0.25, 0.25, 0.25)

NAP  = NAP_calc( CHLz,  0.)
CDOM = CDOM_calc(CHLz, 20.)#10

file_cols = np.vstack((PresCHL, PFT1, PFT2, PFT3, PFT4, CDOM, NAP)).T
np.savetxt(profile_ID + '_IOP.txt', file_cols, header = init_rows, delimiter='\t', comments='')

if len(PresCHL) < 15:
    badstr = '_BAD'
    floatname = profile_ID + '_BAD.nc'
else:
    floatname = profile_ID + '.nc'


################################################################################################  
##########################   phase 4 : Run Fortran code      ###################################
command='./compute.xx ' + profile_ID + '_OASIM.txt ' + profile_ID + '_IOP.txt ' + str(floatname) 
os.system(command)


################################################################################################
#########    phase 5: Prepare irradiance output .nc files for ARGO-model matchup    ############

ncin=NC4.Dataset(floatname,"r")

Ed380_model  =  np.array(ncin.variables['Edz'][3,1:] + ncin.variables['Esz'][3,1:])  * 4 # = 10**(-6) / (10**(-4) * 25) 
Ed412_model  =  np.array(ncin.variables['Edz'][4,1:] + ncin.variables['Esz'][4,1:])  * 4 #  W/m2 to muW/cm2
Ed490_model  =  np.array(ncin.variables['Edz'][7,1:] + ncin.variables['Esz'][7,1:])  * 4

ncin.close()
#Interpolate Ed380 on CHL (OASIM model) depth quotes

Ed380_float = np.interp(PresCHL, Pres380, Ed_380)
Ed412_float = np.interp(PresCHL, Pres412, Ed_412)
Ed490_float = np.interp(PresCHL, Pres490, Ed_490)

ncout = save_matchup(floatname, PresCHL, Ed380_float, Ed412_float, Ed490_float, Ed380_model, Ed412_model, Ed490_model)
################################################################################################


# Copy the in-water radiative transfer model output to a separate directory
movefiles = 'cp ' + str(floatname) + ' NCOUT/'
os.system(movefiles)
