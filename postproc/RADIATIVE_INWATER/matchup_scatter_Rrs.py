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


SIM_MAIN_FOLDER = sys.argv[1]                # SIM_MAIN_FOLDER ='/gpfs/scratch/userexternal/eterzic0/1D_RTM/TESTS/'

CSV_FILE        = open(sys.argv[2], 'r')     # CSV_FILE        = open('../../preproc/RADIATIVE_INWATER/Simulations.csv', 'r')

READER          = csv.reader(CSV_FILE)

fig, ax = plt.subplots(2,3, gridspec_kw = {'wspace':0.25, 'hspace':0.5})
fig.set_size_inches(15,10)

VARLIST           = ['RRS412', 'RRS443', 'RRS490', 'RRS510', 'RRS555' , 'RRS670']

for iline, line in enumerate(READER):  # each line is one simulation
	if iline == 0:
		continue    # we are skipping the header

	SIM_FOLDER = SIM_MAIN_FOLDER + line[0]  # first column of a line

	if not os.path.isdir(SIM_FOLDER + '/RRS'):      # If the MATCHUP folder does not exist, skip it
		continue

	print('Plotting scatter Rrs of simulation %s ... '%line[0], end = '')

	if not os.path.isdir(SIM_MAIN_FOLDER + '/PLOTS'):   
		os.mkdir(SIM_MAIN_FOLDER + '/PLOTS')                # create the PLOTS folder within the simulation folder

	if not os.path.isdir(SIM_MAIN_FOLDER + '/PLOTS/SCATTER'):   
		os.mkdir(SIM_MAIN_FOLDER + '/PLOTS/SCATTER')       # create the SCATTER folder within the PLOTS folder

	if not os.path.isdir(SIM_MAIN_FOLDER + '/PLOTS/SCATTER/Rrs'):   
		os.mkdir(SIM_MAIN_FOLDER + '/PLOTS/SCATTER/Rrs')     # create the Kd folder within the SCATTER folder


	for ifile, filename in enumerate(glob.glob(SIM_FOLDER + '/RRS/*.nc')):

		aux1     = filename.strip('.nc')
		aux2     = aux1.strip(SIM_FOLDER  + '/RRS/')
		wmo, date, lon, lat  = aux2.split('_')
		ncin     = NC4.Dataset(filename, 'r')

		if ifile == 0:
			Rrs_model   = ncin.variables['Rrs'][:,:]
			
		else:
			Rrs_model = np.concatenate((Rrs_model, ncin.variables['Rrs'][:,:]))

		ncin.close()

		Rrs_sat_aux = np.zeros((1,6))

		for ivar, var in enumerate(VARLIST):

			filelist_SAT = SAT_DIR + date.split('-')[0] + "_d-OC_CNR-L3-" + var + "-MedOC4AD4_SAM_1KM-MED-REP-v02.nc"

			ncSAT        = NC4.Dataset(filelist_SAT, 'r')

			lonSAT       = ncSAT.variables['lon'][:]
			latSAT       = ncSAT.variables['lat'][:]

			ilon         = np.argmin(np.abs(lonSAT - lon))
			ilat         = np.argmin(np.abs(latSAT - lat))

			Rrs_sat_aux[0, ivar] = ncSAT.variables[var][0,ilat,ilon].filled(fill_value=np.nan) #np.concatenate((Rrs_sat_aux, [ncSAT.variables[var][0,ilat,ilon].filled(fill_value=np.nan)]))

			ncSAT.close()

		if ifile == 0:
			Rrs_sat     = Rrs_sat_aux
		else:
			Rrs_sat     = np.concatenate((Rrs_sat, Rrs_sat_aux))

	L412 = matchup(Rrs_model[:,0], Rrs_sat[:,0])
	L443 = matchup(Rrs_model[:,1], Rrs_sat[:,1])
	L490 = matchup(Rrs_model[:,2], Rrs_sat[:,2])
	L510 = matchup(Rrs_model[:,3], Rrs_sat[:,3])
	L555 = matchup(Rrs_model[:,4], Rrs_sat[:,4])
	L670 = matchup(Rrs_model[:,5], Rrs_sat[:,5])

	ax1 = plot_matchup_scatter('lin', L412, ax[0,:], 'indigo'  ,     0, 'Rrs', '$Rrs _{\lambda=412}$', 0.60, 0.40, False, None)
	ax2 = plot_matchup_scatter('lin', L443, ax[0,:], 'darkcyan',     1, 'Rrs', '$Rrs _{\lambda=443}$', 0.60, 0.40, False, None)
	ax3 = plot_matchup_scatter('lin', L490, ax[0,:], 'navy'    ,     2, 'Rrs', '$Rrs _{\lambda=490}$', 0.60, 0.40, False, None)
	ax4 = plot_matchup_scatter('lin', L510, ax[1,:], 'forestgreen',  0, 'Rrs', '$Rrs _{\lambda=510}$', 0.60, 0.40, False, None)
	ax5 = plot_matchup_scatter('lin', L555, ax[1,:], 'darkgreen',    1, 'Rrs', '$Rrs _{\lambda=555}$', 0.60, 0.40, False, None)
	ax6 = plot_matchup_scatter('lin', L670, ax[1,:], 'firebrick',    2, 'Rrs', '$Rrs _{\lambda=670}$', 0.60, 0.40, False, None)

	OUTDIR = SIM_MAIN_FOLDER + '/PLOTS/SCATTER/Rrs/' 

	OUTNAME = line[-1].replace(" ", "_").replace(",", "")

	fig.savefig(OUTDIR     + line[0]  + '_plot_' + OUTNAME + '.png', dpi=300)
	print('Saving figure ' + line[0]  + ' plot ' + line[-1])
	sys.stdout.flush()

	ax1.clear()
	ax2.clear()
	ax3.clear()
	ax4.clear()
	ax5.clear()
	ax6.clear()

CSV_FILE.close()
