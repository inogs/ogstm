import glob
from plot_matchup import *
import netCDF4 as NC4

REF_tot = np.zeros((len(glob.glob("*.nc")), )) 

for filename in glob.glob("*.nc"):
#    print(file)
#filename = '6901865_20140830_17.070484_41.44834.nc'

    ncin = NC4.Dataset(filename, 'r')
    ed_float = ncin.variables['Ed_float'][:,:]
    ed_model = ncin.variables['Ed_model'][:,:]
    
    L380 = matchup(ed_model[:,0], ed_float[:,0])
    L412 = matchup(ed_model[:,1], ed_float[:,1])
    L490 = matchup(ed_model[:,2], ed_float[:,2])