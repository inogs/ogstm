# Author Giorgio Bolzon
# useful to check differences between 2 directories
# when md5sum cannot be used because of differences in NetCDF files



import numpy as np
from commons.dataextractor import DataExtractor
from commons.mask import Mask
import os,glob

TheMask=Mask('/gpfs/work/IscrC_MEDCOAST_0/test_swp/TEST01/wrkdir/MODEL/meshmask.nc')
DIR1 = "/gpfs/work/IscrC_MEDCOAST_0/test_swp/test_FSA/ogstm/testcase/TEST01/wrkdir/MODEL/AVE_FREQ_2/"
DIR2 = "/gpfs/work/IscrC_MEDCOAST_0/test_swp/TEST01/wrkdir/MODEL/AVE_FREQ_2/"


searchstr = "ave.20000101-12:00:00*nc"
filelist= glob.glob(DIR1 + searchstr)

for file1 in filelist:
    file2 = DIR2 + os.path.basename(file1) 
    prefix, datestr, varname,_ = os.path.basename(file2).rsplit(".")
    VAR1 = DataExtractor(TheMask,file1,varname).values
    VAR2 = DataExtractor(TheMask,file2,varname).values
    d = VAR2-VAR1
    print varname, (d**2).sum()
    if varname=='O3h': break