import glob
from matchup.statistics import matchup
from plot_matchup import *
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as NC4


runs = ['RUN01/', 'RUN02a/', 'RUN02b/', 'RUN03/', 'RUN04a/', 'RUN04b/', 'RUN04c/', 'RUN04d/', 'RUN04e/']
#pathlist = ['TESTS/RUN01/MATCHUP/',  'TESTS/RUN02b/MATCHUP/', 'TESTS/RUN03/MATCHUP/', 'TESTS/RUN04b/MATCHUP/']
run_name = ['pure_water', 'CDOM min', 'CDOM max','NAP', 'PFT1', 'PFT2', 'PFT3', 'PFT4', 'PFT5']

nWl   = 3
nRuns = len(runs)
nStat = 6

OUTFILE = np.zeros((nWl, nRuns, nStat))

for iRun, run in enumerate(runs):
    
    path = 'TESTS/' + run + 'MATCHUP/'
    
    Ed_float = []
    Ed_model = []

    depth_max = 150.

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
            
    Ed_FLOAT = np.ma.array(Ed_float, mask=(~np.isfinite(Ed_float) | (Ed_float == -999)))
    Ed_MODEL = np.ma.array(Ed_model, mask=(~np.isfinite(Ed_model) | (Ed_model == -999)))
    
    for iWl in range(nWl):

        L = matchup(Ed_MODEL[:,iWl], Ed_FLOAT[:,iWl])
        
        count, bias, RMSE, r_value, slope, intercept = save_stat(L)
        
        OUTFILE[iWl, iRun, 0] = count
        OUTFILE[iWl, iRun, 1] = bias
        OUTFILE[iWl, iRun, 2] = RMSE
        OUTFILE[iWl, iRun, 3] = r_value
        OUTFILE[iWl, iRun, 4] = slope
        OUTFILE[iWl, iRun, 5] = intercept

outdir = '/galileo/home/userexternal/eterzic0/CODE/ogstm/preproc/RADIATIVE_INWATER/STATS/'
np.save('OUTFILE.npy', OUTFILE)
    