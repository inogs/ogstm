from basins import V2 as OGS
from create_IOP import *
from commons.layer import Layer
from commons.time_interval import TimeInterval
from instruments import optbio_float_2019
from instruments import var_conversions
from instruments.matchup_manager import Matchup_Manager
from matchup import statistics
import netCDF4 as NC4
import numpy as np
import os
from profiler import *
import scipy.io.netcdf as NC
from write_matchup import *

command = 'source ../../compilers/machine_modules/galileo.intel'
os.system(command)

wl =np.array( [250., 325., 350., 375., 400.,   425.,  450.,  475.,  500.,  525.,  550.,  575.,  600.,  625.,  650.,  675., 700.,
               725., 775., 850., 950., 1050., 1150., 1250., 1350., 1450., 1550., 1650., 1750., 1900., 2200., 2900., 3700.])
wl_int = wl.astype(np.int64)

wl_str = [np.str(wl_int[i]) for i in range(len(wl_int))] 
str_Ed = ['Ed_' + np.str(wl_str[i]) for i in range(len(wl_str))]
str_Es = ['Es_' + np.str(wl_str[i]) for i in range(len(wl_str))]

maskfile    = '/galileo/home/userexternal/eterzic0/OASIM_POSTPROC/ARGO_MATCHUP/CODES/meshmask.nc'
ncIN=NC.netcdf_file(maskfile,"r")
nav_lev = ncIN.variables['nav_lev'].data.copy()
ncIN.close()

##########################  Create match-up model-data  ##############################
M = Matchup_Manager(ALL_PROFILES,TL,BASEDIR)

TI = TimeInterval("20120101", "20171231","%Y%m%d")
variable='P_l'
varname=var_conversions.FLOAT_OPT_VARS_2019[variable]

Profilelist=optbio_float_2019.FloatSelector(varname,TI , OGS.med)
p = Profilelist[0]
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
CDOM = CDOM_calc(CHLz, 10.)

file_cols = np.vstack((PresCHL, PFT1, PFT2, PFT3, PFT4, CDOM, NAP)).T
np.savetxt(profile_ID + '_IOP.txt', file_cols, header = init_rows, delimiter='\t', comments='')

if len(PresCHL) < 15:
    badstr = '_BAD'
    floatname = profile_ID + '_BAD.nc'
else:
    floatname = profile_ID + '.nc'
    
##########################   phase 4 : Run Fortran code ###########################################
command='./compute.xx ' + profile_ID + '_OASIM.txt ' + profile_ID + '_IOP.txt ' + str(floatname) 
os.system(command)
###################################################################################################

############### phase 5: Prepare irradiance output .nc files for ARGO-model matchup ###############

Ed380_float, Ed380_model, Ed412_float, Ed412_model, Ed490_float, Ed490_model = write_matchup(floatname)

ncout = save_matchup(floatname)

###################################################################################################
# Copy the in-water radiative transfer model output to a separate directory
movefiles = 'cp ' + str(floatname) + ' NCOUT/'
os.system(movefiles)
