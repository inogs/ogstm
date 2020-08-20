#!/bin/env ipython

from __future__ import print_function, division
import datetime
import matplotlib.pyplot as plt
import numpy as np
import csv, glob, os, sys

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

#TI    = TimeInterval('20120101-00:00:00', '20171231-00:00:00', '%Y%m%d-%H:%M:%S')

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

	if not os.path.isdir(SIM_MAIN_FOLDER + '/PLOTS/SCATTER/Kd'):   
		os.mkdir(SIM_MAIN_FOLDER + '/PLOTS/SCATTER/Kd')     # create the Kd folder within the SCATTER folder


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

	TL          = TimeList(datelist)
	TL.filelist = filelist

	for iMonth, month in enumerate(MONTHS):  # Loop over climatological months

        CLIM_MONTH_req       = timerequestors.Clim_month(month)

        ilist, wlist  =  TL.select(CLIM_MONTH_req)

        sum_weight   = 0.
        Kd_model     = np.array([0., 0. , 0.])
        Kd_float     = np.array([0., 0. , 0.])
        std_Kd_model = np.array([0., 0. , 0.])
        std_Kd_float = np.array([0., 0. , 0.])


        for ii, w in zip(ilist,wlist):

			sum_weight += w
			ncin = NC4.Dataset(TL.filelist[ii], 'r')

            Kd_model_aux   = ncin.variables['Kd_model'][:,:]
            Kd_float_aux   = ncin.variables['Kd_float'][:,:]

        	ncin.close()

        	meanval_old   = Kd_model.copy()
			Kd_model     += w/sum_weight*(Kd_model_aux-Kd_model)
			std_Kd_model += w*(Kd_model_aux-meanval_old)*(Kd_model_aux-Kd_model)

        	meanval_old   = Kd_float.copy()
			Kd_float     += w/sum_weight*(Kd_float_aux-Kd_float)
			std_Kd_float += w*(Kd_float_aux-meanval_old)*(Kd_float_aux-Kd_float)


	break













nWl      =  3
nMonths  = 12
nParams  = 10

BOUSSOLE = np.zeros((nWl, nMonths, nParams))
FLOAT    = np.zeros((nWl, nMonths, nParams))

resolution_ls = ['1deg', '2deg', '5deg', 'NWM']

varlist = ['412', '490', 'PAR']

units = [0.01, 0.01, 1.]

unitlist = ['$[\mu W \, m^{-2} \, nm^{-1} ]$', '$[\mu W \, m^{-2} \, nm^{-1} ]$', '$[\mu mol \, quanta \, m^{-2} \, s^{-1}]$']

for resolution in resolution_ls:

	fig, ax = plt.subplots (nrows=3, ncols=1, gridspec_kw = {'wspace':0.25, 'hspace':0.25})
	fig.set_size_inches(9,12)

	MONTHS = np.arange(1, 12 + 1)
	months_str  = ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']

	for iwl, wl in enumerate(varlist):

		float_data = np.load(INPUTDIR_FLOAT + 'OUTFILE_Ed_' + wl.lower() + '_' + resolution + '.npy')

		print(resolution)
		print('OUTFILE_Ed_' + wl.lower() + '_' + resolution + '.npy')

		for im, month in enumerate(MONTHS):
			if wl == 'PAR':
				modelfile = 'ePAR'
			else:
				modelfile = 'ed.' + wl if wl == '490' else 'ed.' + wl + '.5'

			#modelfile += '_%02d.txt'%month
			modelfile += '_%02d_mid_day.txt' % month 

			model_txt = np.loadtxt(INPUTDIR_BOUSSOLE + modelfile, delimiter=' ')

			BOUSSOLE[iwl, im, 0] = model_txt[0] # DATA_mean
			BOUSSOLE[iwl, im, 1] = model_txt[1] # MODEL_mean
			BOUSSOLE[iwl, im, 2] = model_txt[2] # DATA_std
			BOUSSOLE[iwl, im, 3] = model_txt[3] # MODEL_std
			BOUSSOLE[iwl, im, 4] = model_txt[4] # RMSE
			BOUSSOLE[iwl, im, 5] = model_txt[5] # BIAS
			BOUSSOLE[iwl, im, 6] = model_txt[6] # CORR
			BOUSSOLE[iwl, im, 7] = model_txt[7] # SLOPE
			BOUSSOLE[iwl, im, 8] = model_txt[8] # INT
			BOUSSOLE[iwl, im, 9] = model_txt[9] # NO


			FLOAT[iwl, im, 0]    = float_data[im, 0]*units[iwl] # DATA_mean
			FLOAT[iwl, im, 1]    = float_data[im, 1]*units[iwl] # MODEL_mean
			FLOAT[iwl, im, 2]    = float_data[im, 2]*units[iwl] # DATA_std
			FLOAT[iwl, im, 3]    = float_data[im, 3]*units[iwl] # MODEL_std
			FLOAT[iwl, im, 4]    = float_data[im, 4]*units[iwl] # RMSE
			FLOAT[iwl, im, 5]    = float_data[im, 5]*units[iwl] # BIAS
			FLOAT[iwl, im, 6]    = float_data[im, 6]            # CORR
			FLOAT[iwl, im, 7]    = float_data[im, 7]            # SLOPE
			FLOAT[iwl, im, 8]    = float_data[im, 8]            # INT
			FLOAT[iwl, im, 9]    = float_data[im, 9]            # NO


	   
		# Plot mean and standard deviation
		ax[iwl].scatter( MONTHS-0.20, BOUSSOLE[iwl,:,0], s=15,  color='darkblue'      ,     label='B-DATA') 
		ax[iwl].scatter( MONTHS-0.10, BOUSSOLE[iwl,:,1], s=15,  color='dodgerblue'    ,     label='B-MODEL')
		ax[iwl].scatter( MONTHS+0.10, FLOAT[iwl,:,0],    s=15,  color='purple'        ,     label='F-DATA')   
		ax[iwl].scatter( MONTHS+0.20, FLOAT[iwl,:,1],    s=15,  color='palevioletred' ,     label='F-MODEL')
		ax[iwl].errorbar(MONTHS-0.20, BOUSSOLE[iwl,:,0], yerr=BOUSSOLE[iwl,:,2], color='darkblue'      , fmt='o')
		ax[iwl].errorbar(MONTHS-0.10, BOUSSOLE[iwl,:,1], yerr=BOUSSOLE[iwl,:,3], color='dodgerblue'    , fmt='o')
		ax[iwl].errorbar(MONTHS+0.10, FLOAT[iwl,:,0],    yerr=FLOAT[iwl,:,2],    color='purple'        , fmt='o')
		ax[iwl].errorbar(MONTHS+0.20, FLOAT[iwl,:,1],    yerr=FLOAT[iwl,:,3],    color='palevioletred' , fmt='o')
		varname = 'Ed_MEAN_'



		ax[iwl].set_xticks(MONTHS)
		ax[iwl].set_xticklabels(months_str)
		ax[iwl].set_ylabel('Ed ' + wl + unitlist[iwl] )

		ax[iwl].tick_params(axis='both', which='major', labelsize=10)
		
		ax[iwl].set_title('Ed ' + wl, fontsize=16)
	ax[1].legend(loc='lower center', ncol=2, fontsize=12)

	plot_out = '/galileo/home/userexternal/eterzic0/OASIM_POSTPROC/ARGO_MATCHUP_NEW/PLOTS/BOUSSOLE/' + 'Boussole_' + varname + resolution + '_mid_day.png'
	fig.savefig(plot_out, format='png',dpi=150)
	#break




