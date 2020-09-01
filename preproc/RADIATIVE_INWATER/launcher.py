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

CONFIGFILE = open('configuration.txt', 'r')
## Define all the variables you'll read from a .csv file to run the whole set of simulations

fCHL          = float(CONFIGFILE.readline())
fEd           = float(CONFIGFILE.readline())  # fCDOM before
fEs           = float(CONFIGFILE.readline())  # fBBP before
TS_corr       = str2bool(CONFIGFILE.readline().strip('\n')) #bool(CONFIGFILE.readline())
aw_spec       = CONFIGFILE.readline().strip('\n')
a_NAP_model   = CONFIGFILE.readline().strip('\n')
S_NAP         = float(CONFIGFILE.readline())
a_NAP_443     = float(CONFIGFILE.readline())
a_NAP_CORR    = str2bool(CONFIGFILE.readline().strip('\n')) 
a_CDOM_model  = CONFIGFILE.readline().strip('\n')
S_CDOM        = float(CONFIGFILE.readline())
CDOM_TS_corr  = str2bool(CONFIGFILE.readline().strip('\n')) #bool(CONFIGFILE.readline())
aw_380_spec   = CONFIGFILE.readline().strip('\n')
depth_type    = CONFIGFILE.readline().strip('\n')
Kw_type       = CONFIGFILE.readline().strip('\n')
a_CDOM_CORR   = str2bool(CONFIGFILE.readline().strip('\n')) 
a_PFT_use     = str2bool(CONFIGFILE.readline().strip('\n')) #bool(CONFIGFILE.readline())
bp_bbp_model  = CONFIGFILE.readline().strip('\n')
bb_ratio      = float(CONFIGFILE.readline())
bbp_slope     = float(CONFIGFILE.readline())
Q_depth       = CONFIGFILE.readline().strip('\n')
   
CONFIGFILE.close()


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

c_USEFUL = 0

#c_NODATA = 0
#c_LOWED  = 0
c_DEPTH  = 0
#c_CLOUD  = 0
#c_BADED  = 0
c_BADQC  = 0
c_BADKDM = 0
c_BADKDF = 0

for ip in range(ip_start_l,ip_end_l):

	#print("I am %d (%d) running %d (from %d to %d)" % (whoAmI,nWorkers,ip,ip_start_l,ip_end_l))

	# Your serial code goes here
	p = Profilelist[ip]
	profile_ID = p.ID()

	print(profile_ID)

	# Here we don't need the _fitted version of the getMatchups because we are interested only in
	# reading the model ouput's value, which is already at the surface.
	# This is also why we don't need any subset layer!

	List_Ed = [M.getMatchups([p], nav_lev, modelvar, refvar='IRR_380') for modelvar in str_Ed]  # /25 from W m-2 to W m-2 nm -1
	List_Es = [M.getMatchups([p], nav_lev, modelvar, refvar='IRR_380') for modelvar in str_Es]  
	
	# fEd and fEs are multiplication factors for the boundary condition sensitivity analysis
	Ed = np.asarray([0. if len(List_Ed[i].Model)==0 else List_Ed[i].Model[0]/25.*fEd for i in range(len(List_Ed))])  # /25 from W m-2 to W m-2 nm -1
	Es = np.asarray([0. if len(List_Es[i].Model)==0 else List_Es[i].Model[0]/25.*fEs for i in range(len(List_Ed))])  # /25 from W m-2 to W m-2 nm -1
	
	#if np.all(Ed == 0.) and np.all(Es == 0.):
	#	print('I am %d profile %d - No model data for this profile' %(whoAmI, ip))
	#	c_NODATA += 1
	#	continue
	
	#if (Ed[4:9].max() + Es[4:9].max()) < 30./25. :  # it's not 30 but say 30/25 or whatever divided by 25
	#	print('I am %d profile %d - Low irradiance values of OASIM!' %(whoAmI, ip))
	#	c_LOWED += 1
	#	#continue
 
	'''
	phase 2. Read BGC-ARGO profiles
	'''
	PresCHL,   CHLz,    Qc = p.read('CHLA')
	Pres380,   Ed_380,  Qc = p.read('IRR_380')   ; Ed_380 /= 100. # / 100    from mu W cm-2 nm-1 to W m-2 nm -1
	Pres412,   Ed_412,  Qc = p.read('IRR_412')   ; Ed_412 /= 100. # / 100
	Pres490,   Ed_490,  Qc = p.read('IRR_490')   ; Ed_490 /= 100. # / 100
	#PresPAR,   PAR,     Qc = p.read('PAR')
	PresBBP,   BBP700,  Qc = p.read('BBP700')
	PresCDOM,  CDOM,    Qc = p.read('CDOM')
	PresT,     TEMP,    Qc = p.read('TEMP')
	PresS,     SALI,    Qc = p.read('SALI')

	Lon       = p.lon
	Lat       = p.lat
	timestr   = p.time.strftime("%Y%m%d-%H:%M:%S")


	if  np.min([PresCHL[-1], PresBBP[-1], PresCDOM[-1], PresT[-1], PresS[-1]]) < 150. :
		print('I am %d profile %d - Depth range too small' %(whoAmI, ip))
		c_DEPTH += 1
		continue
	
	#if Ed_380[0] < 0.3 or Ed_412[0] < 0.3 or Ed_490[0] < 0.3:
	#	print('I am %d profile %d - BGC-Argo low irradiance values - cloud coverage'  %(whoAmI, ip))
	#	c_CLOUD += 1
	#	continue

	#if np.any(np.diff(Ed_380[0:5]) > 0.) or np.any(np.diff(Ed_412[0:5]) > 0.) or np.any(np.diff(Ed_490[0:5]) > 0.):
	#	print('I am %d profile %d - Increasing RADIOMETRIC values with increasing depth - questionable QC!'  %(whoAmI, ip))
	#	c_BADED += 1
	#	#continue

	''' QC procedures for BGC-Argo profiles '''

	# Create a range vector on which we will interpolate float data for simulations
	Pres = np.arange(0.5, 151.5, 1.)	# You can always change the depth range to 150 m or less - DONE! :-)
	nLevels   = len(Pres)
	init_rows = str(timestr) + '\n' + str(Lat) +  '\n' + str(Lon) + '\n' + str(nLevels)

	CHLz[CHLz < 0.] = 0.  # Check that the CHL profile is placed to zero if negative!!

	# Assign the CHL factor in case you want to run the sensitivity analysis
	CHLz *= fCHL

	CDOM_qc    = CDOM_QC(CDOM)    # before it was multiplied by fCDOM and fBBP for sensitivity analysis purpouses
	BBP700_qc  = BBP700_QC(PresBBP, BBP700)

	# Interpolate to MODEL depth quotes for the simulation runs.
	CHL_int    = np.interp(Pres, PresCHL, CHLz)
	TEMP_int   = np.interp(Pres, PresT, TEMP)
	SALI_int   = np.interp(Pres, PresS, SALI)
	CDOM_int   = np.interp(Pres, PresCDOM, CDOM_qc) 
	BBP700_int = np.interp(Pres, PresBBP, BBP700_qc)
	Ed380_int  = np.interp(Pres, Pres380, Ed_380)
	Ed412_int  = np.interp(Pres, Pres412, Ed_412)
	Ed490_int  = np.interp(Pres, Pres490, Ed_490)

	
	'''
	phase 3. Calculate and save IOPs  
	'''

	####################################################################################################################
	########################################        1 . Pure water         #############################################
	####################################################################################################################    


	if not TS_corr:
		awTS = aw_NO_corr(TEMP_int, SALI_int, model=aw_spec)  
		bwTS = bw_NO_corr(TEMP_int, SALI_int)
	else:
		awTS = aw_TS_corr(TEMP_int, SALI_int, model=aw_spec)  
		bwTS = bw_TS_corr(TEMP_int, SALI_int)                 




	####################################################################################################################
	##################################          2. Non-algal particles - NAP         ################################### 
	####################################################################################################################  
	
	a_NAP = np.zeros((CHL_int.shape[0], wl.shape[0]))

	if a_NAP_model == 'Case1_CHL':
		a_NAP  = aNAP_Case1( CHL_int, CHL_int,  S_NAP)      
	if a_NAP_model == 'Case1_BBP':
		a_NAP  = aNAP_Case1( CHL_int, BBP700_int,  S_NAP)    
	if a_NAP_model == 'Babin_CHL':
		a_NAP = aNAP_Babin(CHL_int, a_NAP_443, S_NAP) 
	if a_NAP_model == 'Babin_BBP':
		a_NAP = aNAP_Babin(BBP700_int, a_NAP_443, S_NAP) 

	if a_NAP_CORR == True:
		a_NAP = aDG_CORR(a_NAP)


	#################################################################################################################### 
	##################################    3. Colored dissolved org. matter - CDOM     ################################## 
	#################################################################################################################### 


	a_CDOM = np.zeros((CHL_int.shape[0], wl.shape[0]))            # If you want to run a simulation without CDOM


	if a_CDOM_model == 'Case1_CHL':
		a_CDOM = aCDOM_Case1(CHL_int,   S_CDOM)     
	if a_CDOM_model == 'Case1_CDOM':
		a_CDOM = aCDOM_Case1_CDOM(CHL_int,  CDOM_int, S_CDOM)
	if a_CDOM_model == 'Kbio_380':

		if not CDOM_TS_corr:
			aw380 = aw_380_NO_corr(TEMP_int, SALI_int, model=aw_380_spec)      # No TS Corr
			bw380 = bw_380_NO_corr(TEMP_int, SALI_int)
		else:
			aw380 = aw_380_TS_corr(TEMP_int, SALI_int, model=aw_380_spec)      # With TS Corr
			bw380 = bw_380_TS_corr(TEMP_int, SALI_int)

		Kbio380, success = calc_Kbio_380(Ed380_int, Pres, TEMP_int, SALI_int, depth_type, aw380, bw380, Kw_type) # MLD or EUPH ; MASON or LIT

		if not success:
			print('I am %d profile %d - BGC-Argo RADIOMETRY QC questionable!!'  %(whoAmI, ip))
			c_BADQC += 1
			continue
		a_CDOM = aCDOM_Kbio(CDOM_int, S_CDOM, Kbio380)

	if a_CDOM_model == 'Kbio_380_CORR':

		if not CDOM_TS_corr:
			aw380 = aw_380_NO_corr(TEMP_int, SALI_int, model=aw_380_spec)      # No TS Corr
			bw380 = bw_380_NO_corr(TEMP_int, SALI_int)
		else:
			aw380 = aw_380_TS_corr(TEMP_int, SALI_int, model=aw_380_spec)      # With TS Corr
			bw380 = bw_380_TS_corr(TEMP_int, SALI_int)

		Kbio380, success = calc_Kbio_380(Ed380_int, Pres, TEMP_int, SALI_int, depth_type, aw380, bw380, Kw_type) # MLD or EUPH ; MASON or LIT

		if not success:
			print('I am %d profile %d - BGC-Argo RADIOMETRY QC questionable!!'  %(whoAmI, ip))
			c_BADQC += 1
			continue

		a_CDOM = aCDOM_Kbio(CDOM_int, S_CDOM, Kbio380)
		a_CDOM[:,:8] -= aNAP_Babin(CDOM_int, a_NAP_443, S_NAP)[:,:8]  # until 500 nm
		if np.any(a_CDOM < 0.): 
			print('I am %d profile %d - NEGATIVE VALUES of aCDOM!!!!' %(whoAmI, ip))

	if a_CDOM_CORR == True:
		a_CDOM = aDG_CORR(a_CDOM)

	#################################################################################################################### 
	##################################    4. Phytoplankton functional types - PFT     ################################## 
	#################################################################################################################### 


	a_PFT_TOT = np.zeros((CHL_int.shape[0], wl.shape[0]))   # If you want to run a simulation without PFTs

	if a_PFT_use == True:
		a_PFT_TOT = PFT_MED(CHL_int)   # That gives us the total aPHY absorption


	#################################################################################################################### 
	###############################    5. Particulate scattering and backscattering   ################################## 
	#################################################################################################################### 

	bp  = np.zeros((CHL_int.shape[0], wl.shape[0]))    # If you want to run a simulation without bp
	bbp = np.zeros((CHL_int.shape[0], wl.shape[0]))

	if bp_bbp_model == 'Case1_CHL':
		bp, bbp  = bp_Case1(CHL_int, bb_ratio)  # backscattering ratio 0.2 to 1.5%  == 0.002 to 0.015
	if bp_bbp_model == 'Case1_BBP':
		bp, bbp  = bp_Case1_bbp(CHL_int, BBP700_int, bb_ratio)   # backscattering ratio  0.2 to 1.5%  == 0.002 to 0.015
	if bp_bbp_model == 'from_BBP700':
		bp,bbp   = bbp_frombbp700(BBP700_int, bbp_slope, bb_ratio)  # slope between 0 and 4. Boss says 1, Organelli uses 2
													  # max bbp is with 0 and 0.015  ; min bbp is with 4 and 0.002


	#################################################################################################################### 
	########################################            SAVE FILES         ############################################# 
	#################################################################################################################### 													  

	write_abw25(wl, awTS, bwTS, fname=profile_ID + '_water_IOP.dat')   

	write_acbc25(wl, a_PFT_TOT, bp, bbp, fname=profile_ID + '_PFT.txt')   # Save PFT_abs, bp and bbp
	
	file_col_DEPTH = Pres.T
	np.savetxt(profile_ID + '_DEPTH.txt', file_col_DEPTH, header = init_rows, delimiter='\t', comments='')
	
	Pres_OUT = Pres.reshape(Pres.shape[0], 1)
	
	file_cols_CDOM = np.hstack((Pres_OUT, a_CDOM))
	np.savetxt(profile_ID + '_CDOM.txt', file_cols_CDOM, delimiter='\t', comments='' )
	
	file_cols_NAP = np.hstack((Pres_OUT, a_NAP))
	np.savetxt(profile_ID + '_NAP.txt',   file_cols_NAP, delimiter='\t', comments='' )
	
	floatname = profile_ID + '.nc'
	
	np.savetxt(profile_ID + '_OASIM.txt', np.c_[Ed, Es])
	
	c_USEFUL += 1

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

	Ed380_model  = INT_wl(375., 400., (ncin.variables['Edz'][3, 1:] + ncin.variables['Esz'][3, 1:]), (ncin.variables['Edz'][4, 1:] + ncin.variables['Esz'][4, 1:]), 380.) #* 4 
	Ed412_model  = INT_wl(400., 425., (ncin.variables['Edz'][4, 1:] + ncin.variables['Esz'][4, 1:]), (ncin.variables['Edz'][5, 1:] + ncin.variables['Esz'][5, 1:]), 412.) #* 4 # = 10**(-6) / (10**(-4) * 25) 
	Ed443_model  = INT_wl(425., 450., (ncin.variables['Edz'][5, 1:] + ncin.variables['Esz'][5, 1:]), (ncin.variables['Edz'][6, 1:] + ncin.variables['Esz'][6, 1:]), 443.) #* 4
	Ed490_model  = INT_wl(475., 500., (ncin.variables['Edz'][7, 1:] + ncin.variables['Esz'][7, 1:]), (ncin.variables['Edz'][8, 1:] + ncin.variables['Esz'][8, 1:]), 490.) #* 4 # W/m2 to muW/cm2 - that was before!
	Ed510_model  = INT_wl(500., 525., (ncin.variables['Edz'][8, 1:] + ncin.variables['Esz'][8, 1:]), (ncin.variables['Edz'][9, 1:] + ncin.variables['Esz'][9, 1:]), 510.) #* 4 # W/m2 to muW/cm2
	Ed555_model  = INT_wl(550., 575., (ncin.variables['Edz'][10,1:] + ncin.variables['Esz'][10,1:]), (ncin.variables['Edz'][11,1:] + ncin.variables['Esz'][11,1:]), 555.) #* 4 # W/m2 to muW/cm2
	Ed670_model  = INT_wl(650., 675., (ncin.variables['Edz'][14,1:] + ncin.variables['Esz'][14,1:]), (ncin.variables['Edz'][15,1:] + ncin.variables['Esz'][15,1:]), 670.) #* 4 # W/m2 to muW/cm2
	

	''' Find the index of the minimum depth range between float and model '''
	
	in380 = np.argmin(np.abs(Pres-Pres380[-1]))
	in412 = np.argmin(np.abs(Pres-Pres412[-1]))
	in490 = np.argmin(np.abs(Pres-Pres490[-1]))
	
	ncout = save_matchup(floatname, Pres[:in380], Pres[:in412], Pres[:in490], Ed380_int[:in380], Ed412_int[:in412], Ed490_int[:in490], Ed380_model[:in380], Ed412_model[:in412], Ed490_model[:in490], timestr)
	
	'''Move the in-water radiative transfer model output to a separate directory'''
	movefiles = 'mv ' + str(floatname) + ' NCOUT/'
	os.system(movefiles)
	

	'''
	phase 6. Calculate Kd from Ed_model and Ed_float at 380, 412 and 490 nm
	'''

	Kd380_model, success380M = calc_Kd(Pres, Ed380_model)
	Kd412_model, success412M = calc_Kd(Pres, Ed412_model)
	Kd490_model, success490M = calc_Kd(Pres, Ed490_model)

	Kd380_float, success380F = calc_Kd(Pres, Ed380_int)
	Kd412_float, success412F = calc_Kd(Pres, Ed412_int)
	Kd490_float, success490F = calc_Kd(Pres, Ed490_int)

	if not success380M or not success412M or not success490M:
		print('I am %d profile %d - Problems computing Kd in MODEL!'  %(whoAmI, ip))
		c_BADKDM += 1
		#continue

	if not success380F or not success412F or not success490F:
		print('I am %d profile %d - Problems computing Kd in FLOAT!'  %(whoAmI, ip))
		c_BADKDF += 1
		#continue

	# Save
	wl_Kd       = [380., 412., 490.]
	Kd_model    = [Kd380_model, Kd412_model, Kd490_model]
	Kd_float    = [Kd380_float, Kd412_float, Kd490_float]

	save_Kd(floatname, Kd_model, Kd_float, wl_Kd, timestr)

	'''
	phase 7. Prepare Ed and Eu values from the model in order to obtain R and then Rrs
	         Satellite wavelenghts for Rrs : 412, 443, 490, 510, 555, 670 - Ed and Eu
	'''

	solz = float(ncin.variables['solz'][0])  # Solar zenih angle in degrees

	Eu412_model  = INT_wl(400., 425., ncin.variables['Euz'][4, 1:], ncin.variables['Euz'][5, 1:], 412.) 
	Eu443_model  = INT_wl(425., 450., ncin.variables['Euz'][5, 1:], ncin.variables['Euz'][6, 1:], 443.) 
	Eu490_model  = INT_wl(475., 500., ncin.variables['Euz'][7, 1:], ncin.variables['Euz'][8, 1:], 490.) 
	Eu510_model  = INT_wl(500., 525., ncin.variables['Euz'][8, 1:], ncin.variables['Euz'][9, 1:], 510.) 
	Eu555_model  = INT_wl(550., 575., ncin.variables['Euz'][10,1:], ncin.variables['Euz'][11,1:], 555.) 
	Eu670_model  = INT_wl(650., 675., ncin.variables['Euz'][14,1:], ncin.variables['Euz'][15,1:], 670.) 
	

	# Compute reflectances

	R412 = Eu412_model / Ed412_model
	R443 = Eu443_model / Ed443_model
	R490 = Eu490_model / Ed490_model
	R510 = Eu510_model / Ed510_model
	R555 = Eu555_model / Ed555_model
	R670 = Eu670_model / Ed670_model


	# Compute Q functions

	Q412 = Q_morel(solz, Pres, CHL_int, 412., Q_depth)
	Q443 = Q_morel(solz, Pres, CHL_int, 443., Q_depth)
	Q490 = Q_morel(solz, Pres, CHL_int, 490., Q_depth)
	Q510 = Q_morel(solz, Pres, CHL_int, 510., Q_depth)
	Q555 = Q_morel(solz, Pres, CHL_int, 555., Q_depth)
	Q670 = Q_morel(solz, Pres, CHL_int, 670., Q_depth)

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


cg_USEFUL = MPI.COMM_WORLD.allreduce(c_USEFUL, op=MPI.SUM)

#cg_NODATA = MPI.COMM_WORLD.allreduce(c_NODATA, op=MPI.SUM)
#cg_LOWED  = MPI.COMM_WORLD.allreduce(c_LOWED , op=MPI.SUM)
cg_DEPTH  = MPI.COMM_WORLD.allreduce(c_DEPTH , op=MPI.SUM)
#cg_CLOUD  = MPI.COMM_WORLD.allreduce(c_CLOUD , op=MPI.SUM)
#cg_BADED  = MPI.COMM_WORLD.allreduce(c_BADED , op=MPI.SUM)
cg_BADQC  = MPI.COMM_WORLD.allreduce(c_BADQC , op=MPI.SUM)
cg_BADKDM  = MPI.COMM_WORLD.allreduce(c_BADKDM , op=MPI.SUM)
cg_BADKDF  = MPI.COMM_WORLD.allreduce(c_BADKDF , op=MPI.SUM)


if whoAmI == 0: 
	print('Number of useful profiles : '                      , cg_USEFUL)
	#print('Number of profiles without model data : '          , cg_NODATA)
	#print('Number of profiles with low model Ed : '           , cg_LOWED)
	print('Number of profiles with low depth range : '        , cg_DEPTH)
	#print('Number of profiles with low float Ed: '            , cg_CLOUD)
	#print('Number of profiles with increasing Ed values: '    , cg_BADED)
	print('Number of profiles with questionable Ed values: '  , cg_BADQC)
	print('Number of profiles with bad Kd MODEL : '           , cg_BADKDM)
	print('Number of profiles with bad Kd FLOAT : '           , cg_BADKDF)




