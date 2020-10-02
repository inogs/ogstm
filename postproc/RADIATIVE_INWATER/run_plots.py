#!/bin/env python

from __future__ import print_function, division
import csv
import glob
import numpy as np
import os, sys

SIM_MAIN_FOLDER = sys.argv[1]                # '/gpfs/scratch/userexternal/eterzic0/1D_RTM/'

for SIM in os.listdir(SIM_MAIN_FOLDER):

	if SIM in ['MAIN', 'BBP']: continue

	FOLDER   = os.path.join(SIM_MAIN_FOLDER, SIM) #+ '/'

	CSV_FILE = glob.glob(FOLDER + '/*.csv')[0]

	print('Plotting simulations %s ... '%SIM, end = '')

	os.system('sbatch plot_Ed.sh %s %s' %(FOLDER, CSV_FILE))
	os.system('sbatch plot_Kd.sh %s %s' %(FOLDER, CSV_FILE))
	os.system('sbatch plot_Rrs.sh %s %s' %(FOLDER, CSV_FILE))
	os.system('sbatch plot_Kd_CLIM.sh %s %s' %(FOLDER, CSV_FILE))
	os.system('sbatch plot_Kd_CLIM_WE.sh %s %s' %(FOLDER, CSV_FILE))
	os.system('sbatch plot_Rrs_CLIM_WE.sh %s %s' %(FOLDER, CSV_FILE))
	os.system('sbatch plot_pcolor_Kd.sh %s %s' %(FOLDER, CSV_FILE))
	os.system('sbatch plot_pcolor_Rrs.sh %s %s' %(FOLDER, CSV_FILE))
	os.system('sbatch plot_pcolor_Kd490_SAT.sh %s %s' %(FOLDER, CSV_FILE))


	print('OK')

