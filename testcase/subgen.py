# generates submask.nc for createGridDA
# to be executed in MODEL/ or where meshmask.nc is

from basins import V0 as OGS
from commons.submask import SubMask
from commons.mask import Mask
TheMask=Mask('meshmask.nc',dzvarname="e3t_0")

for sub in OGS.Pred:
    S=SubMask(sub, TheMask)
    S.save_as_netcdf('submask.nc', maskvarname=sub.name)
