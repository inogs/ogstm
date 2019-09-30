from basins import V2 as OGS
from commons.time_interval import TimeInterval
from instruments import optbio_float_2019
from instruments import var_conversions

TI = TimeInterval("20120101", "20171231","%Y%m%d")
variable='P_l'
varname=var_conversions.FLOAT_OPT_VARS_2019[variable]

Profilelist=optbio_float_2019.FloatSelector(varname,TI , OGS.med)
p = Profilelist[0]

#Read profiles
PresCHL, CHLz,    Qc = p.read('CHL')
Pres380, Ed_380, Qc = p.read('IRR_380')
Pres412, Ed_412, Qc = p.read('IRR_412')
Pres490, Ed_490, Qc = p.read('IRR_490')
PresPAR, PAR,    Qc = p.read('PAR')

Lon = p.lon
Lat = p.lat

timestr = p.time.strftime("%Y%m%d-%H:%M:%S")


command="./edeseu.F90 $CHLz"