#!/bin/env python

from __future__ import print_function, division
import matplotlib.pyplot as plt
import numpy as np
import csv, glob, os, sys

from matchup.statistics import *
import netCDF4 as NC4

from ancillary import *



SIM_MAIN_FOLDER = sys.argv[1]                # '/gpfs/scratch/userexternal/eterzic0/1D_RTM/'

CSV_FILE        = open(sys.argv[2], 'r')     # 'Simulations.csv'

READER          = csv.reader(CSV_FILE)

DEPTH_MAX       = 150.


for iline, line in enumerate(READER):  # each line is one simulation
    if iline == 0:
        continue    # we are skipping the header

    SIM_FOLDER = SIM_MAIN_FOLDER + line[0]  # first column of a line

    if not os.path.isdir(SIM_FOLDER + '/MATCHUP'):      # If the MATCHUP folder does not exist, skip it
        continue

    print('Plotting scatter Ed of simulation %s ... '%line[0], end = '')

    os.mkdir(SIM_FOLDER + '/PLOTS')                # create the PLOTS folder within the simulation folder
    os.mkdir(SIM_FOLDER + '/PLOTS/SCATTER/')       # create the SCATTER folder within the PLOTS folder
    os.mkdir(SIM_FOLDER + '/PLOTS/SCATTER/Ed')     # create the Ed folder within the SCATTER folder

    Ed_float = []
    Ed_model = []

    for filename in glob.glob(SIM_FOLDER + '/MATCHUP/*.nc'):

        aux1 = filename.strip('.nc')
        aux2 = aux1.strip(path)
        wmo, date, lon, lat  = aux2.split('_')
        ncin = NC4.Dataset(filename, 'r')
        depth = ncin.variables['depth'][:]
  
        if depth[0] > 10:
            print('First depth is deeper than 10 m')
            continue
        
        idx = np.where(depth<DEPTH_MAX)

        ed_float = ncin.variables['Ed_float'][idx[0],:]
        ed_model = ncin.variables['Ed_model'][idx[0],:]
        jpk      = ed_float.shape[0]

        for jk in range(jpk):
            Ed_float.append(ed_float[jk,:])
            Ed_model.append(ed_model[jk,:])

    Ed_FLOAT = np.ma.array(Ed_float, mask=(~np.isfinite(Ed_float) | (Ed_float == -999)))
    Ed_MODEL = np.ma.array(Ed_model, mask=(~np.isfinite(Ed_model) | (Ed_model == -999)))

    L380 = matchup(Ed_MODEL[:,0], Ed_FLOAT[:,0])
    L412 = matchup(Ed_MODEL[:,1], Ed_FLOAT[:,1])
    L490 = matchup(Ed_MODEL[:,2], Ed_FLOAT[:,2])

    fig,ax = plt.subplots(1,3, gridspec_kw = {'wspace':0.25, 'hspace':0.5})
    fig.set_size_inches(15,5)

    ax1 = plot_matchup_scatter_Ed('lin', L380, ax, 'indigo'  , 0, '$Ed _{\lambda=380}$', 0.60, 0.40, False, None)
    ax2 = plot_matchup_scatter_Ed('lin', L412, ax, 'darkcyan', 1, '$Ed _{\lambda=412}$', 0.60, 0.40, False, None)
    ax3 = plot_matchup_scatter_Ed('lin', L490, ax, 'navy'    , 2, '$Ed _{\lambda=490}$', 0.60, 0.40, False, None)

    OUTDIR = SIM_FOLDER + '/PLOTS/SCATTER/Ed/' 

    fig.savefig(OUTDIR     + line[0]  + '_plot_' + line[-1].replace(" ", "_") +  '.png', dpi=300)
    print('Saving figure ' + line[0]  + ' plot ' + line[-1])
