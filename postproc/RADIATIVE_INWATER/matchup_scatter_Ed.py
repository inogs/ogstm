#!/bin/env python

from __future__ import print_function, division
import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import csv, glob, os, sys

from matchup.statistics import *
import netCDF4 as NC4

from ancillary import *


SIM_MAIN_FOLDER = sys.argv[1]                # '/gpfs/scratch/userexternal/eterzic0/1D_RTM/'

CSV_FILE        = open(sys.argv[2], 'r')     # 'Simulations.csv'

READER          = csv.reader(CSV_FILE)

fig, ax = plt.subplots(1,3, gridspec_kw = {'wspace':0.25, 'hspace':0.5})

fig.set_size_inches(15,5)


for iline, line in enumerate(READER):  # each line is one simulation
    if iline == 0:
        continue    # we are skipping the header

    SIM_FOLDER = SIM_MAIN_FOLDER + line[0]  # first column of a line

    if not os.path.isdir(SIM_FOLDER + '/MATCHUP'):      # If the MATCHUP folder does not exist, skip it
        continue

    print('Plotting scatter Ed of simulation %s ... '%line[0], end = '')

    if not os.path.isdir(SIM_MAIN_FOLDER + '/PLOTS'):   
        os.mkdir(SIM_MAIN_FOLDER + '/PLOTS')                # create the PLOTS folder within the simulation folder

    if not os.path.isdir(SIM_MAIN_FOLDER + '/PLOTS/SCATTER'):   
        os.mkdir(SIM_MAIN_FOLDER + '/PLOTS/SCATTER')       # create the SCATTER folder within the PLOTS folder

    if not os.path.isdir(SIM_MAIN_FOLDER + '/PLOTS/SCATTER/Ed'):   
        os.mkdir(SIM_MAIN_FOLDER + '/PLOTS/SCATTER/Ed')     # create the Ed folder within the SCATTER folder


    for ifile, filename in enumerate(glob.glob(SIM_FOLDER + '/MATCHUP/*.nc')):


        aux1     = filename.strip('.nc')
        aux2     = aux1.strip(SIM_FOLDER  + '/MATCHUP/')
        wmo, date, lon, lat  = aux2.split('_')
        ncin     = NC4.Dataset(filename, 'r')

        depth380 = ncin.variables['depth380'][:]
        depth412 = ncin.variables['depth412'][:]
        depth490 = ncin.variables['depth490'][:]

        if ifile == 0:

            ed_380   = ncin.variables['Ed_380'][:,:]
            ed_412   = ncin.variables['Ed_412'][:,:]
            ed_490   = ncin.variables['Ed_490'][:,:]

        else:

            ed_380 = np.concatenate((ed_380, ncin.variables['Ed_380'][:,:]))
            ed_412 = np.concatenate((ed_412, ncin.variables['Ed_412'][:,:]))
            ed_490 = np.concatenate((ed_490, ncin.variables['Ed_490'][:,:]))

        ncin.close()

    L380 = matchup(ed_380[:,1], ed_380[:,0])
    L412 = matchup(ed_412[:,1], ed_412[:,0])
    L490 = matchup(ed_490[:,1], ed_490[:,0])

    ax1 = plot_matchup_scatter('lin', L380, ax, 'indigo'  , 0, 'Ed', '$Ed _{\lambda=380}$', 0.60, 0.40, False, None)
    ax2 = plot_matchup_scatter('lin', L412, ax, 'darkcyan', 1, 'Ed', '$Ed _{\lambda=412}$', 0.60, 0.40, False, None)
    ax3 = plot_matchup_scatter('lin', L490, ax, 'navy'    , 2, 'Ed', '$Ed _{\lambda=490}$', 0.60, 0.40, False, None)

    OUTDIR = SIM_MAIN_FOLDER + '/PLOTS/SCATTER/Ed/' 

    OUTNAME = line[-1].replace(" ", "_").replace(",", "")

    fig.savefig(OUTDIR     + line[0]  + '_plot_' + OUTNAME +  '.png', dpi=300)
    print('Saving figure ' + line[0]  + ' plot ' + line[-1])
    sys.stdout.flush()

    ax1.clear()
    ax2.clear()
    ax3.clear()

CSV_FILE.close()
