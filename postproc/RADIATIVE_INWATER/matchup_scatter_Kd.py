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

fig, ax = plt.subplots(1,3, gridspec_kw = {'wspace':0.25, 'hspace':0.5})

fig.set_size_inches(15,5)


for iline, line in enumerate(READER):  # each line is one simulation
    if iline == 0:
        continue    # we are skipping the header

    SIM_FOLDER = SIM_MAIN_FOLDER + line[0]  # first column of a line

    if not os.path.isdir(SIM_FOLDER + '/KD'):      # If the MATCHUP folder does not exist, skip it
        continue

    print('Plotting scatter Ed of simulation %s ... '%line[0], end = '')

    if not os.path.isdir(SIM_MAIN_FOLDER + '/PLOTS'):   
        os.mkdir(SIM_MAIN_FOLDER + '/PLOTS')                # create the PLOTS folder within the simulation folder

    if not os.path.isdir(SIM_MAIN_FOLDER + '/PLOTS/SCATTER'):   
        os.mkdir(SIM_MAIN_FOLDER + '/PLOTS/SCATTER')       # create the SCATTER folder within the PLOTS folder

    if not os.path.isdir(SIM_MAIN_FOLDER + '/PLOTS/SCATTER/Kd'):   
        os.mkdir(SIM_MAIN_FOLDER + '/PLOTS/SCATTER/Kd')     # create the Kd folder within the SCATTER folder


    for ifile, filename in enumerate(glob.glob(SIM_FOLDER + '/KD/*.nc')):


        aux1     = filename.strip('.nc')
        aux2     = aux1.strip(SIM_FOLDER  + '/KD/')
        wmo, date, lon, lat  = aux2.split('_')
        ncin     = NC4.Dataset(filename, 'r')


        if ifile == 0:

            Kd_model   = ncin.variables['Kd_model'][:,:]
            Kd_float   = ncin.variables['Kd_float'][:,:]
            
        else:

            Kd_model = np.concatenate((Kd_model, ncin.variables['Kd_model'][:,:]))
            Kd_float = np.concatenate((Kd_float, ncin.variables['Kd_float'][:,:]))

        ncin.close()

    L380 = matchup(Kd_model[:,0], Kd_float[:,0])
    L412 = matchup(Kd_model[:,1], Kd_float[:,1])
    L490 = matchup(Kd_model[:,2], Kd_float[:,2])

    ax1 = plot_matchup_scatter_Ed('lin', L380, ax, 'indigo'  , 0, '$Kd _{\lambda=380}$', 0.60, 0.40, False, None)
    ax2 = plot_matchup_scatter_Ed('lin', L412, ax, 'darkcyan', 1, '$Kd _{\lambda=412}$', 0.60, 0.40, False, None)
    ax3 = plot_matchup_scatter_Ed('lin', L490, ax, 'navy'    , 2, '$Kd _{\lambda=490}$', 0.60, 0.40, False, None)

    OUTDIR = SIM_MAIN_FOLDER + '/PLOTS/SCATTER/Kd/' 

    OUTNAME = line[-1].replace(" ", "_").replace(",", "")

    fig.savefig(OUTDIR     + line[0]  + '_plot_' + OUTNAME +  '.png', dpi=300)
    print('Saving figure ' + line[0]  + ' plot ' + line[-1])
    sys.stdout.flush()

    ax1.clear()
    ax2.clear()
    ax3.clear()

CSV_FILE.close()
