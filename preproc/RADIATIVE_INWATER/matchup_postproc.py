from basins import V2 as OGS
import glob
from matchup.statistics import matchup
from plot_matchup import *
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as NC4

Ed_FLOAT = []
Ed_MODEL = []

depth_max = 150.

for filename in glob.glob("MATCHUP/*.nc"):
    
    wmo = filename[8:15]
    date = filename[16:24]
    lon = filename[25:34]
    lat = filename[35:43]
    
    ncin = NC4.Dataset(filename, 'r')
    depth = ncin.variables['depth'][:]
    idx = np.where(depth<depth_max)
    ed_float = ncin.variables['Ed_float'][idx[0],:]
    ed_model = ncin.variables['Ed_model'][idx[0],:]
    
    for subbasin in OGS.Pred:
        if subbasin.is_inside(float(lon), float(lat)):
            sub = subbasin
    
    jpk = ed_float.shape[0]
    
    for jk in range(jpk):
        Ed_FLOAT.append(ed_float[jk,:])
        Ed_MODEL.append(ed_model[jk,:])
        
Ed_FLOAT = np.asarray(Ed_FLOAT)
Ed_MODEL = np.asarray(Ed_MODEL)

L380 = matchup(Ed_MODEL[:,0],Ed_FLOAT[:,0])
L412 = matchup(Ed_MODEL[:,1], Ed_FLOAT[:,1])
L490 = matchup(Ed_MODEL[:,2], Ed_FLOAT[:,2])

fig,ax = plt.subplots(1,3, gridspec_kw = {'wspace':0.25, 'hspace':0.5})
fig.set_size_inches(15,5)

ax1 = plot_matchup('log', L380, ax, 'b', 0, 'Ed 380')
ax2 = plot_matchup('log', L412, ax, 'c', 1, 'Ed 412')
ax2 = plot_matchup('log', L490, ax, 'g', 2, 'Ed 490')