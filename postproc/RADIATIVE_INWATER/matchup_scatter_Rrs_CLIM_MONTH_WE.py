#!/bin/env ipython

from __future__ import print_function, division
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


SIM_MAIN_FOLDER = sys.argv[1]                # SIM_MAIN_FOLDER = '/gpfs/scratch/userexternal/eterzic0/1D_RTM/MAIN/'

CSV_FILE        = open(sys.argv[2], 'r')     # CSV_FILE        = open('../../preproc/RADIATIVE_INWATER/Simulations_main.csv', 'r')

READER          = csv.reader(CSV_FILE)

SAT_DIR         = '/gpfs/scratch/userexternal/eterzic0/RRS_OUT/1KM/DAILY/CLIM_MONTH/'

nMonths = 12
MONTHS  = np.arange(1, nMonths + 1)
months_str  = ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']

VARLIST = ['RRS412', 'RRS443', 'RRS490', 'RRS510', 'RRS555' , 'RRS670']


fig1, ax1 = plt.subplots (nrows=3, ncols=1, gridspec_kw = {'wspace':0.25, 'hspace':0.25})
fig1.set_size_inches(12,12)

fig2, ax2 = plt.subplots (nrows=3, ncols=1, gridspec_kw = {'wspace':0.25, 'hspace':0.25})
fig2.set_size_inches(12,12)

for ivar, var in enumerate(VARLIST):

	if ivar == 0:
		Rrs_sat_M_W      = np.load(SAT_DIR + var + '_MEAN_MED_W.npy')
		Rrs_sat_M_E      = np.load(SAT_DIR + var + '_MEAN_MED_E.npy')
		Rrs_sat_S_W      = np.load(SAT_DIR + var + '_STD_MED_W.npy')
		Rrs_sat_S_E      = np.load(SAT_DIR + var + '_STD_MED_E.npy' )
	else:
		Rrs_sat_M_W      = np.vstack((Rrs_sat_M_W, np.load(SAT_DIR + var + '_MEAN_MED_W.npy')))
		Rrs_sat_M_E      = np.vstack((Rrs_sat_M_E, np.load(SAT_DIR + var + '_MEAN_MED_E.npy')))
		Rrs_sat_S_W      = np.vstack((Rrs_sat_S_W, np.load(SAT_DIR + var + '_STD_MED_W.npy' )))
		Rrs_sat_S_E      = np.vstack((Rrs_sat_S_E, np.load(SAT_DIR + var + '_STD_MED_E.npy' )))

for iline, line in enumerate(READER):  # each line is one simulation
	if iline == 0:
		continue    # we are skipping the header

	SIM_FOLDER = SIM_MAIN_FOLDER + line[0]  # first column of a line

	OUTDIR  = SIM_MAIN_FOLDER + '/PLOTS/SCATTER/Rrs_CLIM_MONTH_WE/' 
	OUTNAME = line[-1].replace(" ", "_").replace(",", "")

	if not os.path.isdir(SIM_FOLDER + '/RRS'):      # If the KD folder does not exist, skip it
		continue

	print('Plotting scatter of Rrs - monthly WE of simulation %s ... '%line[0], end = '')

	if not os.path.isdir(SIM_MAIN_FOLDER + '/PLOTS'):   
		os.mkdir(SIM_MAIN_FOLDER + '/PLOTS')                # create the PLOTS folder within the simulation folder

	if not os.path.isdir(SIM_MAIN_FOLDER + '/PLOTS/SCATTER'):   
		os.mkdir(SIM_MAIN_FOLDER + '/PLOTS/SCATTER')        # create the SCATTER folder within the PLOTS folder

	if not os.path.isdir(SIM_MAIN_FOLDER + '/PLOTS/SCATTER/Rrs_CLIM_MONTH_WE'):   
		os.mkdir(SIM_MAIN_FOLDER + '/PLOTS/SCATTER/Rrs_CLIM_MONTH_WE')     # create the Kd folder within the SCATTER folder

	datelist = []
	filelist = []

	for ifile, filename in enumerate(glob.glob(SIM_FOLDER + '/RRS/*.nc')):

		aux1     = filename.strip('.nc')
		aux2     = aux1.strip(SIM_FOLDER  + '/RRS/')
		wmo, date, lon, lat  = aux2.split('_')

		actualtime = datetime.datetime.strptime(date, '%Y%m%d-%H:%M:%S')

		datelist.append(actualtime)
		filelist.append(filename)

	datelist, filelist = (list(t) for t in zip(*sorted(zip(datelist, filelist))))

	TL                 = TimeList(datelist, forceFrequency='daily')
	TL.filelist        = filelist
	TL.timeinterval    = TimeInterval.fromdatetimes(datelist[0], datelist[-1])

	lon                = [float(filename.strip('.nc').strip(SIM_FOLDER + '/RRS/').split('_')[2])  for filename in filelist]
	lat                = [float(filename.strip('.nc').strip(SIM_FOLDER + '/RRS/').split('_')[3])  for filename in filelist]

	for iMonth, month in enumerate(MONTHS):  # Loop over climatological months

		CLIM_MONTH_req = timerequestors.Clim_month(month)

		ilist, wlist   =  TL.select(CLIM_MONTH_req)

		sum_weight_W     = 0.
		sum_weight_E     = 0.

		Rrs_model_W       = np.array([0., 0., 0., 0., 0., 0.])
		Rrs_model_E       = np.array([0., 0., 0., 0., 0., 0.])
		std_Rrs_model_W   = np.array([0., 0., 0., 0., 0., 0.])
		std_Rrs_model_E   = np.array([0., 0., 0., 0., 0., 0.])

		for ii, w in zip(ilist,wlist):

			ncin = NC4.Dataset(TL.filelist[ii], 'r')

			Rrs_model_aux  = ncin.variables['Rrs'][0,:]

			ncin.close()

			if np.any(np.isnan(Rrs_model_aux))  : continue

			if OGS.wes.is_inside(lon[ii], lat[ii]) : 

				sum_weight_W   += w

				meanval_old   = Rrs_model_W.copy()          # value of the contribution to the mean of the previous instant - for the STD
				Rrs_model_W     += w/sum_weight_W*(Rrs_model_aux-meanval_old)                # NewMean = (New-OldMean) * weight/total_weight
				std_Rrs_model_W += w*(Rrs_model_aux-meanval_old)*(Rrs_model_aux-Rrs_model_W)   # NewStd = weight*(New-OldMean)*(New-NewMean) 

			else:

				sum_weight_E   += w

				meanval_old     = Rrs_model_E.copy()          # value of the contribution to the mean of the previous instant - for the STD
				Rrs_model_E     += w/sum_weight_E*(Rrs_model_aux-meanval_old)                # NewMean = (New-OldMean) * weight/total_weight
				std_Rrs_model_E += w*(Rrs_model_aux-meanval_old)*(Rrs_model_aux-Rrs_model_E)   # NewStd = weight*(New-OldMean)*(New-NewMean) 


		std_Rrs_model_W = np.sqrt(std_Rrs_model_W/sum_weight_W)
		std_Rrs_model_E = np.sqrt(std_Rrs_model_E/sum_weight_E)

		if iMonth == 0:

			Rrs_model_M_W  = Rrs_model_W
			Rrs_model_S_W  = std_Rrs_model_W

			Rrs_model_M_E  = Rrs_model_E
			Rrs_model_S_E  = std_Rrs_model_E
			
		else:

			Rrs_model_M_W  = np.vstack((Rrs_model_M_W , Rrs_model_W))
			Rrs_model_S_W  = np.vstack((Rrs_model_S_W , std_Rrs_model_W))

			Rrs_model_M_E  = np.vstack((Rrs_model_M_E , Rrs_model_E))
			Rrs_model_S_E  = np.vstack((Rrs_model_S_E , std_Rrs_model_E))


	for ivar, var in enumerate(VARLIST):

		save_stat(Rrs_model_M_E[:,ivar], Rrs_sat_M_E[ivar,:], OUTDIR + line[0]  + '_stat_' + OUTNAME + '_' + var + '_E.txt' )
		save_stat(Rrs_model_M_W[:,ivar], Rrs_sat_M_W[ivar,:], OUTDIR + line[0]  + '_stat_' + OUTNAME + '_' + var + '_W.txt' )

		# Plot mean and standard deviation

		if ivar < 3:
			ax1[ivar].scatter( MONTHS-0.225, Rrs_model_M_W[:,ivar], s=15,    color='darkblue',        label='MODEL W') 
			ax1[ivar].scatter( MONTHS-0.075, Rrs_sat_M_W[  ivar,:], s=15,    color='dodgerblue',      label='SAT W') 
			ax1[ivar].scatter( MONTHS+0.075, Rrs_model_M_E[:,ivar], s=15,    color='purple'  ,        label='MODEL E') 
			ax1[ivar].scatter( MONTHS+0.225, Rrs_sat_M_E[  ivar,:], s=15,    color='palevioletred',   label='SAT E') 
			
			ax1[ivar].errorbar(MONTHS-0.225, Rrs_model_M_W[:,ivar], yerr=Rrs_model_S_W[:,ivar],  color='darkblue'      , fmt='o')
			ax1[ivar].errorbar(MONTHS-0.075, Rrs_sat_M_W[  ivar,:], yerr=Rrs_sat_S_W[  ivar,:],  color='dodgerblue'    , fmt='o')
			ax1[ivar].errorbar(MONTHS+0.075, Rrs_model_M_E[:,ivar], yerr=Rrs_model_S_E[:,ivar],  color='purple'        , fmt='o')
			ax1[ivar].errorbar(MONTHS+0.225, Rrs_sat_M_E[  ivar,:], yerr=Rrs_sat_S_E[  ivar,:],  color='palevioletred' , fmt='o')

			ax1[ivar].set_xticks(MONTHS)
			ax1[ivar].set_xticklabels(months_str)
			ax1[ivar].set_ylabel(var + ' [$sr^{-1}$]' )
			ax1[ivar].tick_params(axis='both', which='major', labelsize=10)
			ax1[ivar].set_title(var, fontsize=16)
			ax1[ivar].set_ylim(bottom=0)

		else:
			ax2[ivar-3].scatter( MONTHS-0.225, Rrs_model_M_W[:,ivar], s=15,  color='darkblue',      label='MODEL W') 
			ax2[ivar-3].scatter( MONTHS-0.075, Rrs_sat_M_W[  ivar,:], s=15,  color='dodgerblue',    label='SAT W') 
			ax2[ivar-3].scatter( MONTHS+0.075, Rrs_model_M_E[:,ivar], s=15,  color='purple'  ,      label='MODEL E') 
			ax2[ivar-3].scatter( MONTHS+0.225, Rrs_sat_M_E[  ivar,:], s=15,  color='palevioletred', label='SAT E') 

			ax2[ivar-3].errorbar(MONTHS-0.225, Rrs_model_M_W[:,ivar], yerr=Rrs_model_S_W[:,ivar],  color='darkblue'      , fmt='o')
			ax2[ivar-3].errorbar(MONTHS-0.075, Rrs_sat_M_W[  ivar,:], yerr=Rrs_sat_S_W[  ivar,:],  color='dodgerblue'    , fmt='o')
			ax2[ivar-3].errorbar(MONTHS+0.075, Rrs_model_M_E[:,ivar], yerr=Rrs_model_S_E[:,ivar],  color='purple'        , fmt='o')
			ax2[ivar-3].errorbar(MONTHS+0.225, Rrs_sat_M_E[  ivar,:], yerr=Rrs_sat_S_E[  ivar,:],  color='palevioletred' , fmt='o')

			ax2[ivar-3].set_xticks(MONTHS)
			ax2[ivar-3].set_xticklabels(months_str)
			ax2[ivar-3].set_ylabel(var + ' [$sr^{-1}$]' )
			ax2[ivar-3].tick_params(axis='both', which='major', labelsize=10)
			ax2[ivar-3].set_title(var, fontsize=16)
			ax2[ivar-3].set_ylim(bottom=0)

	ax1[2].legend(loc='upper center', ncol=2, fontsize=12)
	ax2[2].legend(loc='upper center', ncol=2, fontsize=12)

	fig1.savefig(OUTDIR     + line[0]  + '_plot_' + OUTNAME +  '_01_WE.png', dpi=300)
	fig2.savefig(OUTDIR     + line[0]  + '_plot_' + OUTNAME +  '_02_WE.png', dpi=300)
	print('Saving figure '  + line[0]  + ' plot ' + line[-1])
	sys.stdout.flush()

	for a in ax1:     
		a.clear()

	for a in ax2:
		a.clear()

CSV_FILE.close()
