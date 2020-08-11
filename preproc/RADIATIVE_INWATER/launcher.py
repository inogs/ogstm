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
from aw_TS_corr import aw_TS_corr, aw_NO_corr, aw_380_TS_corr, aw_380_NO_corr
from bw_TS_corr import bw_TS_corr, bw_NO_corr, bw_380_TS_corr, bw_380_NO_corr
from configuration import *

from commons.layer import Layer
#from instruments.matchup_manager import Matchup_Manager
from matchup_manager_MOD import Matchup_Manager


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

#for ip in range(ip_start_l,ip_end_l):
for ip in range(320,321):

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

	CHLz[CHLz < 0.] = 0.  # Check that the CHL profile is placed to zero if negative!!

	Lon       = p.lon
	Lat       = p.lat
	timestr   = p.time.strftime("%Y%m%d-%H:%M:%S")
	nLevels   = len(PresCHL)
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

	####################################################################################################################
	########################################        1 . Pure water         #############################################
	####################################################################################################################    


	# With T-S correction
	awTS = aw_TS_corr(TEMP_int, SALI_int, model='MASON')  # MASON or LIT
	bwTS = bw_TS_corr(TEMP_int, SALI_int)                 # Then we will change it to a variable

	# Without T-S correction	
	#awTS = aw_NO_corr(TEMP_int, SALI_int, model='MASON')  
	#bwTS = bw_NO_corr(TEMP_int, SALI_int)

	write_abw25(wl, awTS, bwTS, fname=profile_ID + '_water_IOP.dat')   # Save T-S-corrected water IOPs


	####################################################################################################################
	##################################          2. Non-algal particles - NAP         ################################### 
	####################################################################################################################  
	
	
	aNAP  = aNAP_Case1( CHLz,   0.0129)      # 0.0178 max, 0.0104 min and 0.0129 mean


	#################################################################################################################### 
	##################################    3. Colored dissolved org. matter - CDOM     ################################## 
	#################################################################################################################### 

	# 3.1. Case 1 type
	aCDOM = aCDOM_Case1(np.zeros(CHLz.shape),   0.017)     # If you want to run a simulation without CDOM
	#aCDOM = aCDOM_Case1(CHLz,   0.017)     # 0.02   max, 0.015  min and 0.017  mean 

	CDOM_qc = CDOM_QC(CDOM) # quality-controlled CDOM profile
	CDOM_int = np.interp(PresCHL, PresCDOM, CDOM_qc) # Interpolate CDOM to CHL depth!
	# You need this because at the moment you're saving all IOPs on CHLz depth quotas for the model run.

	# 3.2. Case 1 type with CDOM shape
	#aCDOM = aCDOM_Case1_CDOM(CHLz,  CDOM_int, 0.017)     # 0.02   max, 0.015  min and 0.017  mean 

	# 3.3. Estimating aCDOM(380) from Kd380

	#aw380 = aw_380_NO_corr(TEMP_int, SALI_int, model='LIT')    # No TS Corr
	#bw380 = bw_380_NO_corr(TEMP_int, SALI_int)

	#aw380 = aw_380_TS_corr(TEMP_int, SALI_int, model='MASON')    # With TS Corr
	#bw380 = bw_380_TS_corr(TEMP_int, SALI_int)

	#Kbio380 = calc_Kbio_380(Ed_380, Pres380, PresCHL, TEMP_int, SALI_int, 'EUPH', aw380, bw380, 'MASON') # MLD or EUPH ; MASON or LIT

	#aCDOM = aCDOM_Kbio(CDOM_int, 0.017, Kbio380)

	#################################################################################################################### 
	##################################    4. Phytoplankton functional types - PFT     ################################## 
	#################################################################################################################### 


	aPFT_TOT = np.zeros((CHLz.shape[0], wl.shape[0]))   # If you want to run a simulation without PFTs
	#aPFT_TOT = PFT_MED(CHLz)   # That gives us the total aPHY absorption


	#################################################################################################################### 
	###############################    5. Particulate scattering and backscattering   ################################## 
	#################################################################################################################### 


	BBP700_qc  = BBP700_QC(PresBBP, BBP700)
	BBP700_int = np.interp(PresCHL, PresBBP, BBP700_qc) # Interpolate bbp700 to CHL depth!

	bp  = np.zeros((CHLz.shape[0], wl.shape[0]))    # If you want to run a simulation without bp
	bbp = np.zeros((CHLz.shape[0], wl.shape[0]))

	# 5.1. Case 1 - Chl shape
	#bp, bbp  = bp_Case1(CHLz, 0.002)  # backscattering ratio 0.2 to 1.5%  == 0.002 to 0.015

	# 5.2 Case 1 - bbp shape
	#bp, bbp  = bp_Case1_bbp(CHLz, BBP700_int, 0.002)   # backscattering ratio  0.2 to 1.5%  == 0.002 to 0.015

	# 5.3 bp, bbp from bbp700
	#bp,bbp   = bbp_frombbp700(BBP700_int, 4., 0.002)  # slope between 0 and 4. Boss says 1, Organelli uses 2
													  # max bbp is with 0 and 0.015  ; min bbp is with 4 and 0.002


	#################################################################################################################### 
	########################################            SAVE FILES         ############################################# 
	#################################################################################################################### 													  


	write_acbc25(wl, aPFT_TOT, bp, bbp, fname=profile_ID + '_PFT.txt')   # Save PFT_abs, bp and bbp
	
	file_col_DEPTH = PresCHL.T
	np.savetxt(profile_ID + '_DEPTH.txt', file_col_DEPTH, header = init_rows, delimiter='\t', comments='')
	
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
	command='./compute.xx ' + profile_ID + '_OASIM.txt ' + profile_ID + '_water_IOP.dat ' + profile_ID + '_DEPTH.txt ' + profile_ID + '_PFT.txt '  + profile_ID + '_CDOM.txt '  + profile_ID + '_NAP.txt ' + str(floatname) + ' >> log'

	print ('I am %d profile %d - %s ' %(whoAmI, ip,command ))
	subprocess.call(command, shell=True)
	
	'''
	phase 5: Prepare irradiance output .nc files for ARGO-model matchup
	'''  
	ncin=NC4.Dataset(floatname,"r")

	Ed380_model  = INT_wl(375., 400., (ncin.variables['Edz'][3, 1:] + ncin.variables['Esz'][3, 1:]), (ncin.variables['Edz'][4, 1:] + ncin.variables['Esz'][4, 1:]), 380.) * 4 
	Ed412_model  = INT_wl(400., 425., (ncin.variables['Edz'][4, 1:] + ncin.variables['Esz'][4, 1:]), (ncin.variables['Edz'][5, 1:] + ncin.variables['Esz'][5, 1:]), 412.) * 4 # = 10**(-6) / (10**(-4) * 25) 
	Ed443_model  = INT_wl(425., 450., (ncin.variables['Edz'][5, 1:] + ncin.variables['Esz'][5, 1:]), (ncin.variables['Edz'][6, 1:] + ncin.variables['Esz'][6, 1:]), 443.) * 4
	Ed490_model  = INT_wl(475., 500., (ncin.variables['Edz'][7, 1:] + ncin.variables['Esz'][7, 1:]), (ncin.variables['Edz'][8, 1:] + ncin.variables['Esz'][8, 1:]), 490.) * 4 # W/m2 to muW/cm2
	Ed510_model  = INT_wl(500., 525., (ncin.variables['Edz'][8, 1:] + ncin.variables['Esz'][8, 1:]), (ncin.variables['Edz'][9, 1:] + ncin.variables['Esz'][9, 1:]), 510.) * 4 # W/m2 to muW/cm2
	Ed555_model  = INT_wl(550., 575., (ncin.variables['Edz'][10,1:] + ncin.variables['Esz'][10,1:]), (ncin.variables['Edz'][11,1:] + ncin.variables['Esz'][11,1:]), 555.) * 4 # W/m2 to muW/cm2
	Ed670_model  = INT_wl(650., 675., (ncin.variables['Edz'][14,1:] + ncin.variables['Esz'][14,1:]), (ncin.variables['Edz'][15,1:] + ncin.variables['Esz'][15,1:]), 670.) * 4 # W/m2 to muW/cm2
	

	'''Interpolate Ed380 on CHL (OASIM model) depth quotes'''
	
	Ed380_float = np.interp(PresCHL, Pres380, Ed_380)
	Ed412_float = np.interp(PresCHL, Pres412, Ed_412)
	Ed490_float = np.interp(PresCHL, Pres490, Ed_490)
	
	ncout = save_matchup(floatname, PresCHL, Ed380_float, Ed412_float, Ed490_float, Ed380_model, Ed412_model, Ed490_model, timestr)
	
	'''Move the in-water radiative transfer model output to a separate directory'''
	movefiles = 'mv ' + str(floatname) + ' NCOUT/'
	os.system(movefiles)
	

	'''
	phase 6. Prepare Ed and Eu values from the model in order to obtain R and then Rrs
	         Satellite wavelenghts for Rrs : 412, 443, 490, 510, 555, 670 - Ed and Eu
	'''

	solz = float(ncin.variables['solz'][0])  # Solar zenih angle in degrees

	Eu412_model  = INT_wl(400., 425., ncin.variables['Euz'][4, 1:], ncin.variables['Euz'][5, 1:], 412.) * 4 # = 10**(-6) / (10**(-4) * 25) 
	Eu443_model  = INT_wl(425., 450., ncin.variables['Euz'][5, 1:], ncin.variables['Euz'][6, 1:], 443.) * 4
	Eu490_model  = INT_wl(475., 500., ncin.variables['Euz'][7, 1:], ncin.variables['Euz'][8, 1:], 490.) * 4 # W/m2 to muW/cm2
	Eu510_model  = INT_wl(500., 525., ncin.variables['Euz'][8, 1:], ncin.variables['Euz'][9, 1:], 510.) * 4 # W/m2 to muW/cm2
	Eu555_model  = INT_wl(550., 575., ncin.variables['Euz'][10,1:], ncin.variables['Euz'][11,1:], 555.) * 4 # W/m2 to muW/cm2
	Eu670_model  = INT_wl(650., 675., ncin.variables['Euz'][14,1:], ncin.variables['Euz'][15,1:], 670.) * 4 # W/m2 to muW/cm2
	

	# Compute reflectances

	R412 = Eu412_model / Ed412_model
	R443 = Eu443_model / Ed443_model
	R490 = Eu490_model / Ed490_model
	R510 = Eu510_model / Ed510_model
	R555 = Eu555_model / Ed555_model
	R670 = Eu670_model / Ed670_model


	# Compute Q functions

	Q412 = Q_morel(solz, CHLz, 412.)
	Q443 = Q_morel(solz, CHLz, 443.)
	Q490 = Q_morel(solz, CHLz, 490.)
	Q510 = Q_morel(solz, CHLz, 510.)
	Q555 = Q_morel(solz, CHLz, 555.)
	Q670 = Q_morel(solz, CHLz, 670.)

	Q = [Q412, Q443, Q490, Q510, Q555, Q670]

	print('I am %d profile %d - Q values' %(whoAmI, ip) , Q)

	# Remote sensing reflectances

	Rrs412 = Q412 * R412[0]
	Rrs443 = Q443 * R443[0]
	Rrs490 = Q490 * R490[0]
	Rrs510 = Q510 * R510[0]
	Rrs555 = Q555 * R555[0]
	Rrs670 = Q670 * R670[0]

	wl_RRS = [412., 443., 490., 510., 555., 670.]
	Rrs    = [Rrs412, Rrs443, Rrs490, Rrs510, Rrs555, Rrs670]

	ncRrs = save_reflectance(floatname, wl_RRS, Rrs, timestr)
	
	''' Move the .txt files you don't need any more '''
	txtfiles1 = 'mv ' + profile_ID + '_OASIM.txt'     + ' TXT_FILES/' 
	txtfiles2 = 'mv ' + profile_ID + '_PFT.txt'       + ' TXT_FILES/' 
	txtfiles3 = 'mv ' + profile_ID + '_NAP.txt'       + ' TXT_FILES/' 
	txtfiles4 = 'mv ' + profile_ID + '_CDOM.txt'      + ' TXT_FILES/'
	txtfiles5 = 'mv ' + profile_ID + '_water_IOP.dat' + ' TXT_FILES/'
	txtfiles6 = 'mv ' + profile_ID + '_DEPTH.txt'     + ' TXT_FILES/'
	
	ncin.close()

	os.system(txtfiles1)
	os.system(txtfiles2)
	os.system(txtfiles3)
	os.system(txtfiles4)
	os.system(txtfiles5)
	os.system(txtfiles6)

