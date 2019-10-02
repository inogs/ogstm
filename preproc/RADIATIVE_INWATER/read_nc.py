from basins import V2 as OGS
from commons.layer import Layer
from commons.time_interval import TimeInterval
from instruments import optbio_float_2019
from instruments import var_conversions
import os
from profiler import *
import scipy.io.netcdf as NC

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

L_380 = M.getMatchups([p], nav_lev, 'Ed_380').subset(Layer(0,1.5))
L_412 = M.getMatchups([p], nav_lev, 'Ed_412').subset(Layer(0,1.5))
L_490 = M.getMatchups([p], nav_lev, 'Ed_490').subset(Layer(0,1.5))


#Read profiles
PresCHL, CHLz,    Qc = p.read('CHL')
Pres380, Ed_380, Qc = p.read('IRR_380')
Pres412, Ed_412, Qc = p.read('IRR_412')
Pres490, Ed_490, Qc = p.read('IRR_490')
PresPAR, PAR,    Qc = p.read('PAR')

Lon = p.lon
Lat = p.lat

timestr = p.time.strftime("%Y%m%d-%H:%M:%S")

#pres_CHL, CHL_pft1, CHL_pft2, CHL_pft3, CHL_pft4 CDOM NAP IOP.txt

# .txt OASIM -matchup

#command="./compute.xx " + timestr + ' ' + Lat 
#os.system(command)

