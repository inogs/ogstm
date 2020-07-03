from basins import V2 as OGS
import glob
from matchup.statistics import matchup
from plot_matchup import *
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as NC4

Ed_FLOAT_W = []   ; Ed_FLOAT_E = []  
Ed_MODEL_W = []   ; Ed_MODEL_E = [] 

depth_max = 150.

for filename in glob.glob("MATCHUP/*.nc"):
    
    aux1 = filename.strip('.nc')
    aux2 = aux1.strip('MATCHUP/')
    
    wmo, date, lon, lat  = aux2.split('_')
    
    ncin = NC4.Dataset(filename, 'r')
    depth = ncin.variables['depth'][:]
    
    if depth[0] > 10:
        continue
    idx = np.where(depth<depth_max)
    ed_float = ncin.variables['Ed_float'][idx[0],:]
    ed_model = ncin.variables['Ed_model'][idx[0],:]
    
    jpk = ed_float.shape[0]
    
    for subbasin in OGS.wes: # eas2
        if subbasin.is_inside(float(lon), float(lat)):
            for jk in range(jpk):
                Ed_FLOAT_W.append(ed_float[jk,:])
                Ed_MODEL_W.append(ed_model[jk,:])
        else:
            for jk in range(jpk):
                Ed_FLOAT_E.append(ed_float[jk,:])
                Ed_MODEL_E.append(ed_model[jk,:])
        
Ed_FLOAT_W = np.asarray(Ed_FLOAT_W)   ; Ed_FLOAT_E = np.asarray(Ed_FLOAT_E)

Ed_MODEL_W = np.asarray(Ed_MODEL_W)   ; Ed_MODEL_E = np.asarray(Ed_MODEL_E)

L380_W = matchup(Ed_MODEL_W[:,0], Ed_FLOAT_W[:,0]) ; L380_E = matchup(Ed_MODEL_E[:,0], Ed_FLOAT_E[:,0])
L412_W = matchup(Ed_MODEL_W[:,1], Ed_FLOAT_W[:,1]) ; L412_E = matchup(Ed_MODEL_E[:,1], Ed_FLOAT_E[:,1])
L490_W = matchup(Ed_MODEL_W[:,2], Ed_FLOAT_W[:,2]) ; L490_E = matchup(Ed_MODEL_E[:,2], Ed_FLOAT_E[:,2])

fig,ax = plt.subplots(1,3, gridspec_kw = {'wspace':0.25, 'hspace':0.5})
fig.set_size_inches(15,5)


ax1 = plot_matchup('lin', L380_W, ax, 'indigo', 0, 'Ed 380', 0.70, 0.30, True, 'West')
ax1 = plot_matchup('lin', L380_E, ax, 'mediumpurple', 0, 'Ed 380', 0.10, 0.70, True, 'East')

ax2 = plot_matchup('lin', L412_W, ax, 'darkcyan', 1, 'Ed 412', 0.70, 0.30, True, 'West')
ax2 = plot_matchup('lin', L412_E, ax, 'mediumaquamarine', 1, 'Ed 412', 0.10, 0.70, True, 'East')

ax3 = plot_matchup('lin', L490_W, ax, 'navy', 2, 'Ed 490', 0.70, 0.30, True, 'West')
ax3 = plot_matchup('lin', L490_E, ax, 'lightsteelblue', 2, 'Ed 490', 0.10, 0.70, True, 'East')
