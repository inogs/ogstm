import glob
from matchup.statistics import matchup
from plot_matchup import *
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as NC4

Ed_float = []
Ed_model = []

depth_max = 150.
path = 'TESTS/RUN01/MATCHUP/'

for filename in glob.glob(path + "*.nc"):
    
    aux1 = filename.strip('.nc')
    aux2 = aux1.strip(path)
    
    wmo, date, lon, lat  = aux2.split('_')
    
    ncin = NC4.Dataset(filename, 'r')
    depth = ncin.variables['depth'][:]
    
    if depth[0] > 10:
        continue
    idx = np.where(depth<depth_max)
    ed_float = ncin.variables['Ed_float'][idx[0],:]
    ed_model = ncin.variables['Ed_model'][idx[0],:]
    
    jpk = ed_float.shape[0]
    
    for jk in range(jpk):
        Ed_float.append(ed_float[jk,:])
        Ed_model.append(ed_model[jk,:])
        
Ed_FLOAT = np.ma.array(Ed_float, mask=(~np.isfinite(Ed_float) | (Ed_float == -999)))#np.asarray(Ed_FLOAT)
Ed_MODEL = np.ma.array(Ed_model, mask=(~np.isfinite(Ed_model) | (Ed_model == -999)))#np.asarray(Ed_MODEL)

L380 = matchup(Ed_MODEL[:,0],Ed_FLOAT[:,0])
L412 = matchup(Ed_MODEL[:,1], Ed_FLOAT[:,1])
L490 = matchup(Ed_MODEL[:,2], Ed_FLOAT[:,2])

fig,ax = plt.subplots(1,3, gridspec_kw = {'wspace':0.25, 'hspace':0.5})
fig.set_size_inches(15,5)

ax1 = plot_matchup('lin', L380, ax, 'indigo', 0, 'Ed 380', 0.70, 0.30,False, None)
ax2 = plot_matchup('lin', L412, ax, 'darkcyan', 1, 'Ed 412', 0.70, 0.30, False, None)
ax3 = plot_matchup('lin', L490, ax, 'navy', 2, 'Ed 490', 0.70, 0.30, False, None)