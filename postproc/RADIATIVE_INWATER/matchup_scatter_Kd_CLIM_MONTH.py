#!/bin/env ipython

from __future__ import print_function, division
import datetime
import matplotlib.pyplot as plt
import numpy as np
import csv, glob, os, sys

from commons import timerequestors
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from matchup.statistics import *
import netCDF4 as NC4

from ancillary import *


SIM_MAIN_FOLDER = sys.argv[1]                # '/gpfs/scratch/userexternal/eterzic0/1D_RTM/TESTS/'

CSV_FILE        = open(sys.argv[2], 'r')     # 'Simulations.csv'

READER          = csv.reader(CSV_FILE)

nMonths = 12
MONTHS  = np.arange(1, nMonths + 1)

VARLIST = ['380', '412', '490']

fig, ax = plt.subplots (nrows=3, ncols=1, gridspec_kw = {'wspace':0.25, 'hspace':0.25})
fig.set_size_inches(9,12)


for iline, line in enumerate(READER):  # each line is one simulation
	if iline == 0:
		continue    # we are skipping the header

	SIM_FOLDER = SIM_MAIN_FOLDER + line[0]  # first column of a line

	if not os.path.isdir(SIM_FOLDER + '/KD'):      # If the KD folder does not exist, skip it
		continue

	print('Plotting scatter Kd of simulation %s ... '%line[0], end = '')

	if not os.path.isdir(SIM_MAIN_FOLDER + '/PLOTS'):   
		os.mkdir(SIM_MAIN_FOLDER + '/PLOTS')                # create the PLOTS folder within the simulation folder

	if not os.path.isdir(SIM_MAIN_FOLDER + '/PLOTS/SCATTER'):   
		os.mkdir(SIM_MAIN_FOLDER + '/PLOTS/SCATTER')        # create the SCATTER folder within the PLOTS folder

	if not os.path.isdir(SIM_MAIN_FOLDER + '/PLOTS/SCATTER/Kd_CLIM_MONTH'):   
		os.mkdir(SIM_MAIN_FOLDER + '/PLOTS/SCATTER/Kd_CLIM_MONTH')     # create the Kd folder within the SCATTER folder


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

	TL                 = TimeList(datelist)
	TL.filelist        = filelist

	for iMonth, month in enumerate(MONTHS):  # Loop over climatological months

		CLIM_MONTH_req = timerequestors.Clim_month(month)

		ilist, wlist   =  TL.select(CLIM_MONTH_req)

		sum_weight     = 0.
		Kd_model       = np.array([0., 0., 0.])
		Kd_float       = np.array([0., 0., 0.])
		std_Kd_model   = np.array([0., 0., 0.])
		std_Kd_float   = np.array([0., 0., 0.])

		for ii, w in zip(ilist,wlist):

			sum_weight += w
			ncin = NC4.Dataset(TL.filelist[ii], 'r')

			Kd_model_aux  = ncin.variables['Kd_model'][0,:]
			Kd_float_aux  = ncin.variables['Kd_float'][0,:]

			ncin.close()

			meanval_old   = Kd_model.copy()          # value of the contribution to the mean of the previous instant - for the STD
			Kd_model     += w/sum_weight*(Kd_model_aux-meanval_old)                # NewMean = (New-OldMean) * weight/total_weight
			std_Kd_model += w*(Kd_model_aux-meanval_old)*(Kd_model_aux-Kd_model)   # NewStd = weight*(New-OldMean)*(New-NewMean) 

			meanval_old   = Kd_float.copy()
			Kd_float     += w/sum_weight*(Kd_float_aux-meanval_old)
			std_Kd_float += w*(Kd_float_aux-meanval_old)*(Kd_float_aux-Kd_float)

		std_Kd_model = np.sqrt(std_Kd_model/sum_weight)
		std_Kd_float = np.sqrt(std_Kd_float/sum_weight)
	
	for ivar, var in enumerate(VARLIST):
	   
		# Plot mean and standard deviation
		ax[ivar].scatter( MONTHS-0.15, Kd_model[ivar], s=15,                    color='darkblue', label='MODEL') 
		ax[ivar].scatter( MONTHS+0.15, Kd_float[ivar], s=15,                    color='purple'  , label='FLOAT')   
		ax[ivar].errorbar(MONTHS-0.15, Kd_model[ivar], yerr=std_Kd_model[ivar], color='darkblue', fmt='o')
		ax[ivar].errorbar(MONTHS+0.15, Kd_float[ivar], yerr=std_Kd_float[ivar], color='purple'  , fmt='o')
		varname = 'Kd_MEAN_'

		ax[ivar].set_xticks(MONTHS)
		ax[ivar].set_xticklabels(months_str)
		ax[ivar].set_ylabel('Kd ' + var + unitlist[ivar] )

		ax[ivar].tick_params(axis='both', which='major', labelsize=10)
		
		ax[ivar].set_title('Kd ' + var, fontsize=16)

	ax[1].legend(loc='lower center', ncol=2, fontsize=12)


	OUTDIR  = SIM_MAIN_FOLDER + '/PLOTS/SCATTER/Kd_CLIM_MONTH/' 
	OUTNAME = line[-1].replace(" ", "_").replace(",", "")
	fig.savefig(OUTDIR     + line[0]  + '_plot_' + OUTNAME +  '.png', dpi=300)
	print('Saving figure ' + line[0]  + ' plot ' + line[-1])
	sys.stdout.flush()

	for a in ax:     
		a.clear()

CSV_FILE.close()
