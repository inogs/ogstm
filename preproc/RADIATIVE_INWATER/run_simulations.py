#!/bin/env python

from __future__ import print_function, division
import csv
import numpy as np
import os, sys

SIM_MAIN_FOLDER = '/gpfs/scratch/userexternal/eterzic0/1D_RTM/'

CSV_FILE = open('Simulations.csv', 'r')     #sys.argv[1] 

READER   = csv.reader(CSV_FILE)

for iline, line in enumerate(READER):  # each line is one simulation
	if iline == 0:
		continue    # we are skipping the header

	SIM_FOLDER = SIM_MAIN_FOLDER + line[0]  # first column of a line

	if os.path.isdir(SIM_FOLDER):      # Skip the directory if there is already a simulation inside
		continue

	os.mkdir(SIM_FOLDER)               # create the simulation folder
	os.mkdir(SIM_FOLDER + '/CODE')     # create the CODE folder within the simulation folder

	os.system('cp *.py %s/CODE'%SIM_FOLDER)
	os.system('cp compute.xx %s/CODE'%SIM_FOLDER)
	os.system('cp Profilelist.pkl %s/CODE'%SIM_FOLDER)
	os.system('cp runit2.sh %s/CODE'%SIM_FOLDER)
	os.system('cp -pr bcs %s/CODE'%SIM_FOLDER)
	os.system('cp env.sh %s/CODE'%SIM_FOLDER)
	os.system('cp Sullivan_T_chart.txt %s/CODE'%SIM_FOLDER)


	CONFIG_FILE = open(SIM_FOLDER + '/CODE/configuration.txt')

	CONFIG_FILE.write('%f\n'%float(line[1]))   # fCHL
	CONFIG_FILE.write('%f\n'%float(line[2]))   # fCDOM
	CONFIG_FILE.write('%f\n'%float(line[3]))   # fBBP
	CONFIG_FILE.write('%s\n'%(line[4]))        # TS_corr
	CONFIG_FILE.write('%s\n'%(line[5]))        # aw_spec
	CONFIG_FILE.write('%s\n'%(line[6]))        # a_NAP_model
	CONFIG_FILE.write('%f\n'%float(line[7]))   # a_NAP_443
	CONFIG_FILE.write('%f\n'%float(line[8]))   # S_NAP
	CONFIG_FILE.write('%s\n'%(line[9]))        # a_CDOM_model
	CONFIG_FILE.write('%f\n'%float(line[10]))  # S_CDOM
	CONFIG_FILE.write('%s\n'%(line[11]))       # CDOM_TS_corr
	CONFIG_FILE.write('%s\n'%(line[12]))       # aw_380_spec
	CONFIG_FILE.write('%s\n'%(line[13]))       # depth_type
	CONFIG_FILE.write('%s\n'%(line[14]))       # Kw_type
	CONFIG_FILE.write('%s\n'%(line[15]))       # a_PFT_use
	CONFIG_FILE.write('%s\n'%(line[16]))       # bp_bbp_model
	CONFIG_FILE.write('%f\n'%float(line[17]))  # bb_ratio
	CONFIG_FILE.write('%f\n'%float(line[18]))  # bbp_slope
	CONFIG_FILE.write('%s\n'%(line[19]))       # Q_depth

	CONFIG_FILE.close()

	break