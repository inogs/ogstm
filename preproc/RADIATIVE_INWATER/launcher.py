#!/bin/env python
#
# Launcher using python MPI4Py
#
# It expects the user to provide either two inputs 
#
from __future__ import print_function, division
import sys, numpy as np
from mpi4py import MPI
import os
import pprint, pickle
import subprocess

from ancillary import *
from aw_TS_corr import aw_TS_corr, aw_NO_corr
from bw_TS_corr import bw_TS_corr, bw_NO_corr
from configuration import *

from commons.layer import Layer
from instruments.matchup_manager import Matchup_Manager


## MPI Ancillary functions ##
comm     = MPI.COMM_WORLD  # Communications macro
whoAmI   = comm.Get_rank() # Who are you? who? who?
nWorkers = comm.Get_size() # Total number of processors used (workers)

def worksplit(istart,iend):
	'''
	Divide the work between the processors
	'''
	istart_l, iend_l = istart, iend
	irange = iend - istart
	if (nWorkers < irange):
		# We split normally among processes assuming no remainder
		rangePerProcess = int(np.floor(irange/nWorkers))
		istart_l = istart   + whoAmI*rangePerProcess
		iend_l   = istart_l + rangePerProcess
		# Handle the remainder
		remainder = irange - rangePerProcess*nWorkers
		if remainder > whoAmI:
			istart_l += whoAmI
			iend_l   += whoAmI+1;
		else:
			istart_l += remainder
			iend_l   += remainder
	else:
		# Each process will forcefully conduct one instant.
		istart_l = whoAmI   if whoAmI < iend else iend
		iend_l   = whoAmI+1 if whoAmI < iend else iend

	return istart_l, iend_l	


## Main script ##
# Here is where you unpickle your profile list
pkl_file = open('Profilelist.pkl', 'rb')

Profilelist = pickle.load(pkl_file)
Floatlist   = pickle.load(pkl_file)

pkl_file.close()


MAXPROFILES = len(Profilelist) # = 100 
ip_start, ip_end = -1, -1

# Parse command line arguments
if len(sys.argv) == 1:
	# No arguments provided, hence take the maximum of the list
	ip_start, ip_end = 0, MAXPROFILES
if len(sys.argv) == 3:
	# User has provided start and end points
	ip_start, ip_end = int(sys.argv[1]), int(sys.argv[2])
	if ip_start < 0 or ip_start > MAXPROFILES: raise ValueError("Wrong start point %d (%d)" % (ip_start,MAXPROFILES))
	if   ip_end < 0 or   ip_end > MAXPROFILES: raise ValueError("Wrong end   point %d (%d)" % (ip_end,MAXPROFILES))

# Wrong inputs
if ip_start < 0 or ip_end < 0: raise ValueError("Wrong number of input arguments!")

# Now each rank should have the list of profiles and where the user
# wants to start and end globally (ip_start, p_end). 
# Each processor (worker) will loop a subset of that list (ip_start_l, ip_end_l).
# These are determined by the worksplit
ip_start_l, ip_end_l = worksplit(ip_start,ip_end)

# Each processor (worker) loops from ip_start_l to ip_end_l (+1 for python reasons...)
func = lambda Pres, E0, k : E0 * np.exp(-k*Pres)

M = Matchup_Manager(Profilelist,TL,BASEDIR)

for ip in range(ip_start_l,ip_end_l):
#for ip in range(320,321):

	#print("I am %d (%d) running %d (from %d to %d)" % (whoAmI,nWorkers,ip,ip_start_l,ip_end_l))
	# Your serial code goes here

	p = Profilelist[ip]
	profile_ID = p.ID()

	print(profile_ID)
	
	List_Ed = [M.getMatchups_fitted([p], nav_lev, modelvar, func, refvar='IRR_380').subset(Layer(0,0.1)) for modelvar in str_Ed]
	List_Es = [M.getMatchups_fitted([p], nav_lev, modelvar, func, refvar='IRR_380').subset(Layer(0,0.1)) for modelvar in str_Es]
	
	Ed = np.asarray([0. if len(List_Ed[i].Model)==0 else List_Ed[i].Model[0] for i in range(len(List_Ed))])
	Es = np.asarray([0. if len(List_Es[i].Model)==0 else List_Es[i].Model[0] for i in range(len(List_Ed))])
	
	if Ed.all() == 0. and Es.all() == 0.:
		print('I am %d profile %d - No model data for this profile' %(whoAmI, ip))
		continue
	
	if (Ed[4:9].max() + Es[4:9].max()) < 30.:
		print('I am %d profile %d - Low irradiance values of OASIM!' %(whoAmI, ip))
		continue
 
	'''
	phase 2. Read BGC-ARGO profiles
	'''
	PresCHL,   CHLz,    Qc = p.read('CHLA')
	Pres380,   Ed_380,  Qc = p.read('IRR_380')
	Pres412,   Ed_412,  Qc = p.read('IRR_412')
	Pres490,   Ed_490,  Qc = p.read('IRR_490')
	#PresPAR,   PAR,     Qc = p.read('PAR')
	PresBBP,   BBP700,  Qc = p.read('BBP700')
	PresCDOM,  CDOM,    Qc = p.read('CDOM')
	PresT,     TEMP,    Qc = p.read('TEMP')
	PresS,     SALI,    Qc = p.read('SALI')
	#PresN,     NIT,     Qc = p.read('NITRATE')
	#PresO,     DOXY,    Qc = p.read('DOXY')


	CHLz[CHLz < 0.] = 0.

	Lon = p.lon
	Lat = p.lat
	timestr = p.time.strftime("%Y%m%d-%H:%M:%S")
	nLevels = len(PresCHL)
	init_rows = str(timestr) + '\n' + str(Lat) +  '\n' + str(Lon) + '\n' + str(nLevels)
	
	if PresCHL[0] == 0. :#or PresCDOM[0] == 0. or PresBBP[0] == 0.:
		print('I am %d profile %d - First depth equals 0' %(whoAmI, ip))
		continue
	
	if Ed_380[0] < 30. or Ed_412[0] < 30. or Ed_490[0] < 30.:
		print('I am %d profile %d - BGC-Argo low irradiance values - cloud coverage'  %(whoAmI, ip))
		continue
	
	if PresCHL.max() < np.min([Pres380.max(), Pres412.max(), Pres490.max()]):
		print('I am %d profile %d - Depth range too small' %(whoAmI, ip))
		continue

	if PresCHL.shape[0] < np.max([PresT.shape[0], PresS.shape[0]]):
		print('I am %d profile %d - Cannot interpolate to CHL depth quotes' %(whoAmI, ip))
		continue

	# Interpolate TEMP and SALI to CHL depth quotes
	TEMP_int = np.interp(PresCHL, PresT, TEMP)
	SALI_int = np.interp(PresCHL, PresS, SALI)
	
	'''
	phase 3. Calculate and save IOPs  
	'''

	# Pure water

	# With T-S correction
	awTS = aw_TS_corr(TEMP_int, SALI_int, model='MASON')  # MASON or LIT
	bwTS = bw_TS_corr(TEMP_int, SALI_int)                 # Then we will change it to a variable

	#awTS = aw_NO_corr(TEMP_int, SALI_int, model='MASON')  # Without T-S correction
	#bwTS = bw_NO_corr(TEMP_int, SALI_int)

	write_abw25(wl, awTS, bwTS, fname=profile_ID + '_water_IOP.dat')   # Save T-S-corrected water IOPs

	PFT1, PFT2, PFT3, PFT4 = PFT_calc(CHLz, 0., 0., 0., 0.)  #0.40, 0.30, 0.25, 0.05)
	
	aNAP  = aNAP_Case1( CHLz,   0.0129) #0.0104 min , 0.0178 max and 0.0129 mean
	aCDOM = aCDOM_Case1(CHLz,   0.)
	
	file_cols_PFT = np.vstack((PresCHL, PFT1, PFT2, PFT3, PFT4)).T
	np.savetxt(profile_ID + '_PFT.txt', file_cols_PFT, header = init_rows, delimiter='\t', comments='')
	
	Pres = PresCHL.reshape(PresCHL.shape[0], 1)
	
	file_cols_CDOM = np.hstack((Pres, aCDOM))
	np.savetxt(profile_ID + '_CDOM.txt', file_cols_CDOM, delimiter='\t', comments='' )
	
	file_cols_NAP = np.hstack((Pres, aNAP))
	np.savetxt(profile_ID + '_NAP.txt',   file_cols_NAP, delimiter='\t', comments='' )
	
	floatname = profile_ID + '.nc'
	
	np.savetxt(profile_ID + '_OASIM.txt', np.c_[Ed, Es])
	
	
	'''  
	phase 4 : Run Fortran code
	'''
	command='./compute.xx ' + profile_ID + '_OASIM.txt ' + profile_ID + '_water_IOP.dat ' + profile_ID + '_PFT.txt '  + profile_ID + '_CDOM.txt '  + profile_ID + '_NAP.txt ' + str(floatname) + ' >> log'

	print ('I am %d profile %d - %s ' %(whoAmI, ip,command ))
	subprocess.call(command, shell=True)
	
	'''
	phase 5: Prepare irradiance output .nc files for ARGO-model matchup
	'''  
	ncin=NC4.Dataset(floatname,"r")
	
	Ed380_model  =  np.array( 0.8  * (ncin.variables['Edz'][3,1:] + ncin.variables['Esz'][3,1:])  + 0.2 *  (ncin.variables['Edz'][4,1:] + ncin.variables['Esz'][4,1:]))  * 4 # = 10**(-6) / (10**(-4) * 25) 
	Ed412_model  =  np.array( 0.52 * (ncin.variables['Edz'][4,1:] + ncin.variables['Esz'][4,1:])  + 0.48 * (ncin.variables['Edz'][5,1:] + ncin.variables['Esz'][5,1:]))  * 4 #  W/m2 to muW/cm2
	Ed490_model  =  np.array( 0.4  * (ncin.variables['Edz'][7,1:] + ncin.variables['Esz'][7,1:])  + 0.6 *  (ncin.variables['Edz'][8,1:] + ncin.variables['Esz'][8,1:]))  * 4
	
	ncin.close()
	'''Interpolate Ed380 on CHL (OASIM model) depth quotes'''
	
	Ed380_float = np.interp(PresCHL, Pres380, Ed_380)
	Ed412_float = np.interp(PresCHL, Pres412, Ed_412)
	Ed490_float = np.interp(PresCHL, Pres490, Ed_490)
	
	ncout = save_matchup(floatname, PresCHL, Ed380_float, Ed412_float, Ed490_float, Ed380_model, Ed412_model, Ed490_model, timestr)
	
	'''Move the in-water radiative transfer model output to a separate directory'''
	movefiles = 'mv ' + str(floatname) + ' NCOUT/'
	os.system(movefiles)
	
	
	''' Move the .txt files you don't need any more '''
	txtfiles1 = 'mv ' + profile_ID + '_OASIM.txt'     + ' TXT_FILES/' 
	txtfiles2 = 'mv ' + profile_ID + '_PFT.txt'       + ' TXT_FILES/' 
	txtfiles3 = 'mv ' + profile_ID + '_NAP.txt'       + ' TXT_FILES/' 
	txtfiles4 = 'mv ' + profile_ID + '_CDOM.txt'      + ' TXT_FILES/'
	txtfiles5 = 'mv ' + profile_ID + '_water_IOP.dat' + ' TXT_FILES/'
	os.system(txtfiles1)
	os.system(txtfiles2)
	os.system(txtfiles3)
	os.system(txtfiles4)
	os.system(txtfiles5)
