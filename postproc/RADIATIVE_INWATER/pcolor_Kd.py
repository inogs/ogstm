#!/bin/env ipython
from __future__ import print_function
import datetime
import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import csv, glob, os, sys

from basins import V2 as OGS
from commons import timerequestors
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from matchup.statistics import *
import netCDF4 as NC4

from ancillary import *


SIM_MAIN_FOLDER = sys.argv[1]                # SIM_MAIN_FOLDER = '/gpfs/scratch/userexternal/eterzic0/1D_RTM/TESTS/'

CSV_FILE        = open(sys.argv[2], 'r')     # CSV_FILE        = open('../../preproc/RADIATIVE_INWATER/Simulations.csv', 'r')

READER          = csv.reader(CSV_FILE)

#SAT_DIR         = '/gpfs/scratch/userexternal/eterzic0/KD490_OUT/1KM/DAILY/CLIM_MONTH/'

#Kd_sat_M_W      = np.load(SAT_DIR + 'Kd_MEAN_MED_W.npy')
#Kd_sat_M_E      = np.load(SAT_DIR + 'Kd_MEAN_MED_E.npy')
#Kd_sat_S_W      = np.load(SAT_DIR + 'Kd_STD_MED_W.npy')
#Kd_sat_S_E      = np.load(SAT_DIR + 'Kd_STD_MED_E.npy')

nMonths     = 12
MONTHS      = np.arange(1, nMonths + 1)
months_str  = ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']
nSub        = len(OGS.med.basin_list)
nStat       = 4

VARLIST = ['380', '412', '490']

basin_list_abbrev = ['ALB', 'SWM1', 'SWM2', 'NWM', 'TYR1', 'TYR2', 'ADR1', 'ADR2', 'AEG', 'ION1', 'ION2', 'ION3', 'LEV1', 'LEV2', 'LEV3', 'LEV4' ]

fig1,ax1 = plt.subplots(2,2)
fig2,ax2 = plt.subplots(2,2)
fig3,ax3 = plt.subplots(2,2)

fig1.set_size_inches(11.69, 8.27)
fig2.set_size_inches(11.69, 8.27)
fig3.set_size_inches(11.69, 8.27)


for iline, line in enumerate(READER):  # each line is one simulation
    if iline == 0:
        continue    # we are skipping the header

    SIM_FOLDER = SIM_MAIN_FOLDER + line[0]  # first column of a line

    if not os.path.isdir(SIM_FOLDER + '/KD'):      # If the KD folder does not exist, skip it
        continue

    print('Plotting scatter Kd of simulation %s ... '%line[0], end = '')

    if not os.path.isdir(SIM_MAIN_FOLDER + '/PLOTS'):   
        os.mkdir(SIM_MAIN_FOLDER + '/PLOTS')                # create the PLOTS folder within the simulation folder

    if not os.path.isdir(SIM_MAIN_FOLDER + '/PLOTS/PCOLOR'):   
        os.mkdir(SIM_MAIN_FOLDER + '/PLOTS/PCOLOR')        # create the SCATTER folder within the PLOTS folder

    if not os.path.isdir(SIM_MAIN_FOLDER + '/PLOTS/PCOLOR/Kd_CLIM'):   
        os.mkdir(SIM_MAIN_FOLDER + '/PLOTS/PCOLOR/Kd_CLIM')     # create the Kd folder within the SCATTER folder

    datelist = []
    filelist = []

    for ifile, filename in enumerate(glob.glob(SIM_FOLDER + '/KD/*.nc')):

        aux1     = filename.strip('.nc')
        aux2     = aux1.strip(SIM_FOLDER  + '/KD/')
        wmo, date, lon, lat  = aux2.split('_')

        actualtime = datetime.datetime.strptime(date, '%Y%m%d-%H:%M:%S')

        datelist.append(actualtime)
        filelist.append(filename)

    datelist, filelist = (list(t) for t in zip(*sorted(zip(datelist, filelist))))

    TL                 = TimeList(datelist, forceFrequency='daily')
    TL.filelist        = filelist
    TL.timeinterval    = TimeInterval.fromdatetimes(datelist[0], datelist[-1])

    lon                = [float(filename.strip('.nc').strip(SIM_FOLDER + '/KD/').split('_')[2])  for filename in filelist]
    lat                = [float(filename.strip('.nc').strip(SIM_FOLDER + '/KD/').split('_')[3])  for filename in filelist]

    Kd_380_SUB = np.zeros((nSub, nMonths, nStat))
    Kd_412_SUB = np.zeros((nSub, nMonths, nStat))
    Kd_490_SUB = np.zeros((nSub, nMonths, nStat))


    for isub, sub in enumerate(OGS.med):         # Loop over subbasins

        for iMonth, month in enumerate(MONTHS):  # Loop over climatological months

            CLIM_MONTH_req = timerequestors.Clim_month(month)

            ilist, wlist   =  TL.select(CLIM_MONTH_req)

            for ii, w in zip(ilist,wlist):

                if not sub.is_inside(lon[ii], lat[ii]): continue

                ncin = NC4.Dataset(TL.filelist[ii], 'r')

                if ii == 0:
                    Kd_model   = ncin.variables['Kd_model'][:,:]
                    Kd_float   = ncin.variables['Kd_float'][:,:]
    
                else:
                    Kd_model = np.concatenate((Kd_model, ncin.variables['Kd_model'][:,:]))
                    Kd_float = np.concatenate((Kd_float, ncin.variables['Kd_float'][:,:]))

                ncin.close()

            L380 = matchup(Kd_model[:,0], Kd_float[:,0])
            L412 = matchup(Kd_model[:,1], Kd_float[:,1])
            L490 = matchup(Kd_model[:,2], Kd_float[:,2])


            Kd_380_SUB[isub, iMonth, :] = [L380.Model.mean(), L380.Ref.mean(), L380.bias(), L380.RMSE() ]
            Kd_412_SUB[isub, iMonth, :] = [L412.Model.mean(), L412.Ref.mean(), L412.bias(), L412.RMSE() ]
            Kd_490_SUB[isub, iMonth, :] = [L490.Model.mean(), L490.Ref.mean(), L490.bias(), L490.RMSE() ]


    plot_pcolor(ax1, Kd_380_SUB[:,:,0], Kd_380_SUB[:,:,1], Kd_380_SUB[:,:,2], Kd_380_SUB[:,:,3], '380', basin_list_abbrev, months_str)
    plot_pcolor(ax2, Kd_412_SUB[:,:,0], Kd_412_SUB[:,:,1], Kd_412_SUB[:,:,2], Kd_412_SUB[:,:,3], '412', basin_list_abbrev, months_str)
    plot_pcolor(ax3, Kd_490_SUB[:,:,0], Kd_490_SUB[:,:,1], Kd_490_SUB[:,:,2], Kd_490_SUB[:,:,3], '490', basin_list_abbrev, months_str)

    fig1.tight_layout()
    fig2.tight_layout()
    fig3.tight_layout()

    OUTDIR = SIM_MAIN_FOLDER + '/PLOTS/PCOLOR/Kd_CLIM/' 

    OUTNAME = line[-1].replace(" ", "_").replace(",", "")

    fig1.savefig(OUTDIR     + line[0]  + '_plot_380_' + OUTNAME +  '.png', dpi=300)
    fig2.savefig(OUTDIR     + line[0]  + '_plot_412_' + OUTNAME +  '.png', dpi=300)
    fig3.savefig(OUTDIR     + line[0]  + '_plot_490_' + OUTNAME +  '.png', dpi=300)

    print('Saving figure ' + line[0]  + ' plot ' + line[-1]
    sys.stdout.flush()

    ax1.clear()
    ax2.clear()
    ax3.clear()

CSV_FILE.close()

    
