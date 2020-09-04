#!/bin/env ipython

from __future__ import print_function, division
import datetime
import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import csv, glob, os, sys

from commons import timerequestors
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from matchup.statistics import *
import netCDF4 as NC4

from ancillary import *


SIM_MAIN_FOLDER = sys.argv[1]                # SIM_MAIN_FOLDER = '/gpfs/scratch/userexternal/eterzic0/1D_RTM/TESTS/'

CSV_FILE        = open(sys.argv[2], 'r')     # CSV_FILE        = open('../../preproc/RADIATIVE_INWATER/Simulations.csv', 'r')

READER          = csv.reader(CSV_FILE)

SAT_DIR         = '/gpfs/scratch/userexternal/eterzic0/KD490_OUT/1KM/DAILY/CLIM_MONTH/'

Kd_sat_M        = np.load(SAT_DIR + 'Kd_MEAN_MED.npy')
Kd_sat_S        = np.load(SAT_DIR + 'Kd_STD_MED.npy')

nMonths = 12
MONTHS  = np.arange(1, nMonths + 1)
months_str  = ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']

VARLIST = ['380', '412', '490']

fig, ax = plt.subplots (nrows=3, ncols=1, gridspec_kw = {'wspace':0.25, 'hspace':0.25})
fig.set_size_inches(9,12)


for iline, line in enumerate(READER):  # each line is one simulation
	if iline == 0:
		continue    # we are skipping the header

	SIM_FOLDER = SIM_MAIN_FOLDER + line[0]  # first column of a line

	if not os.path.isdir(SIM_FOLDER + '/KD'):      # If the KD folder does not exist, skip it
		continue

	print('Plotting scatter Kd CLIM MONTH of simulation %s ... '%line[0], end = '')

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


	TL                 = TimeList(datelist, forceFrequency='daily')
	TL.filelist        = filelist
	TL.timeinterval    = TimeInterval.fromdatetimes(datelist[0], datelist[-1])


	for iMonth, month in enumerate(MONTHS):  # Loop over climatological months

		CLIM_MONTH_req = timerequestors.Clim_month(month)

		ilist, wlist   =  TL.select(CLIM_MONTH_req)

		sum_weight     = 0.
		Kd_model       = np.array([0., 0., 0.])
		Kd_float       = np.array([0., 0., 0.])
		std_Kd_model   = np.array([0., 0., 0.])
		std_Kd_float   = np.array([0., 0., 0.])

		for ii, w in zip(ilist,wlist):

			ncin = NC4.Dataset(TL.filelist[ii], 'r')

			Kd_model_aux  = ncin.variables['Kd_model'][0,:]
			Kd_float_aux  = ncin.variables['Kd_float'][0,:]

			ncin.close()

			if np.any(np.isnan(Kd_model_aux)) or np.any(np.isnan(Kd_float_aux)) : continue

			sum_weight   += w

			meanval_old   = Kd_model.copy()          # value of the contribution to the mean of the previous instant - for the STD
			Kd_model     += w/sum_weight*(Kd_model_aux-meanval_old)                # NewMean = (New-OldMean) * weight/total_weight
			std_Kd_model += w*(Kd_model_aux-meanval_old)*(Kd_model_aux-Kd_model)   # NewStd = weight*(New-OldMean)*(New-NewMean) 

			meanval_old   = Kd_float.copy()
			Kd_float     += w/sum_weight*(Kd_float_aux-meanval_old)
			std_Kd_float += w*(Kd_float_aux-meanval_old)*(Kd_float_aux-Kd_float)

		std_Kd_model = np.sqrt(std_Kd_model/sum_weight)
		std_Kd_float = np.sqrt(std_Kd_float/sum_weight)

		if iMonth == 0:

			Kd_model_M  = Kd_model
			Kd_float_M  = Kd_float
			Kd_model_S  = std_Kd_model
			Kd_float_S  = std_Kd_float
			
		else:

			Kd_model_M  = np.vstack((Kd_model_M , Kd_model))
			Kd_float_M  = np.vstack((Kd_float_M , Kd_float))
			Kd_model_S  = np.vstack((Kd_model_S , std_Kd_model))
			Kd_float_S  = np.vstack((Kd_float_S , std_Kd_float))

	
	for ivar, var in enumerate(VARLIST):
	   
	   	if ivar < 2:
		# Plot mean and standard deviation
			ax[ivar].scatter( MONTHS-0.15, Kd_model_M[:,ivar], s=15,                    color='dodgerblue', label='MODEL')  # or darkblue
			ax[ivar].scatter( MONTHS+0.15, Kd_float_M[:,ivar], s=15,                    color='purple'    , label='FLOAT')    # or palevioletred
			
			ax[ivar].errorbar(MONTHS-0.15, Kd_model_M[:,ivar], yerr=Kd_model_S[:,ivar], color='dodgerblue', fmt='o')
			ax[ivar].errorbar(MONTHS+0.15, Kd_float_M[:,ivar], yerr=Kd_float_S[:,ivar], color='purple'    , fmt='o')

		else:
			ax[ivar].scatter( MONTHS-0.20, Kd_model_M[:,ivar], s=15,                    color='dodgerblue', label='MODEL') # or darkblue
			ax[ivar].scatter( MONTHS     , Kd_sat_M[:],        s=15,                    color='peru' , label='SAT')    # or palevioletred
			ax[ivar].scatter( MONTHS+0.20, Kd_float_M[:,ivar], s=15,                    color='purple'    , label='FLOAT')    # or palevioletred
			
			ax[ivar].errorbar(MONTHS-0.20, Kd_model_M[:,ivar], yerr=Kd_model_S[:,ivar], color='dodgerblue', fmt='o')
			ax[ivar].errorbar(MONTHS,      Kd_sat_M[:],        yerr=Kd_sat_S[:],        color='peru' , fmt='o') # or sandybrown, before limegreen
			ax[ivar].errorbar(MONTHS+0.20, Kd_float_M[:,ivar], yerr=Kd_float_S[:,ivar], color='purple'    , fmt='o')


		ax[ivar].set_xticks(MONTHS)
		ax[ivar].set_xticklabels(months_str)
		ax[ivar].set_ylabel('Kd ' + var + ' [$m^{-1}$]' )

		ax[ivar].tick_params(axis='both', which='major', labelsize=10)
		
		ax[ivar].set_title('Kd ' + var, fontsize=16)

	ax[2].legend(loc='upper center', ncol=2, fontsize=12)


	OUTDIR  = SIM_MAIN_FOLDER + '/PLOTS/SCATTER/Kd_CLIM_MONTH/' 
	OUTNAME = line[-1].replace(" ", "_").replace(",", "")
	fig.savefig(OUTDIR     + line[0]  + '_plot_' + OUTNAME +  '.png', dpi=300)
	print('Saving figure ' + line[0]  + ' plot ' + line[-1])
	sys.stdout.flush()

	for a in ax:     
		a.clear()

CSV_FILE.close()
