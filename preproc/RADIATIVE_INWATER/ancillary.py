#!/bin/env python

import netCDF4 as NC4
import numpy as np
from scipy import ndimage

from configuration import wl

''' Here you put all the functions, also for the IOPs , T-S corrections, etc. '''

def write_abw25(wl, aw, bw, fname='bcs/abw25.dat'):

	file  = open(fname, 'w')
	for idepth in range(aw.shape[1]):
		for iwl in range(aw.shape[0]):
			file.write('%5d%15.4f%10.4f\n'%(wl[iwl], aw[iwl, idepth], bw[iwl, idepth]))
	file.close()

	return 

def running_mean(a,WSZ):
		'''
		Smoothes a 1-D numpy array.
		
		WSZ: smoothing window size needs, which must be odd number,
		as in the original MATLAB implementation.
		'''
		import numpy as np
		out0  = np.convolve(a,np.ones(WSZ,dtype=int),'valid')/WSZ    
		r     = np.arange(1,WSZ-1,2)
		start = np.cumsum(a[:WSZ-1])[::2]/r
		stop  = (np.cumsum(a[:-WSZ:-1])[::2]/r)[::-1]
		return np.concatenate((  start , out0, stop  ))


def findVars(Varlist, allvars=['CHLA', 'IRR_380', 'IRR_412', 'IRR_490', 'TEMP']):
	if len(Varlist)==0: return False    
	
	for var in allvars:
		if not var in Varlist:
			return False
	return True

'''
Prepare model and BGC-Argo output for .nc files and match-up analysis
'''

def save_matchup(ncfile, PresCHL, Ed380_float, Ed412_float, Ed490_float, Ed380_model, Ed412_model, Ed490_model, timestr):

	modelfile = 'MATCHUP/' + ncfile
	ncmodel   = NC4.Dataset(modelfile,"w");
				
	ncdepth = ncmodel.createDimension('depth',     len(PresCHL));
	ncwave  = ncmodel.createDimension('wavelength', 3);

	setattr(ncmodel, 'time', timestr);
	
	ncDepth = ncmodel.createVariable('depth', 'f', ('depth')); 
	setattr(ncDepth, 'unit',  '[m]' );
	ncDepth[:] = PresCHL
	
	ncEdf = ncmodel.createVariable('Ed_float', 'f', ('depth', 'wavelength'));
	setattr(ncEdf, 'missing_value',-1.0 );     
	setattr(ncEdf, 'long_name',  'Downward irradiance ' );     
	setattr(ncEdf, 'unit',  '[uW/cm2/nm]' );

	ncEdf[:] = np.vstack((Ed380_float, Ed412_float, Ed490_float)).T
	
	ncEdm = ncmodel.createVariable('Ed_model', 'f', ('depth', 'wavelength'));
	setattr(ncEdm, 'missing_value',-1.0 );     
	setattr(ncEdm, 'long_name',  'Downward irradiance ' );     
	setattr(ncEdm, 'unit',  '[uW/cm2/nm]' );

	ncEdm[:] = np.vstack((Ed380_model, Ed412_model, Ed490_model)).T
	
	ncmodel.close()

	return ncmodel

# BGC-Argo profile QC procedures

def CDOM_QC(CDOM):
	CDOM_MF = ndimage.median_filter(CDOM, size=5)
	CDOM_QC = running_mean(CDOM_QC, 7)
	return CDOM_QC

def BBP700_QC(PresBBP, BBP700):
	BBP700_MF = ndimage.median_filter(BBP700, size=11)
	BBP700_RM = running_mean(BBP700_MF, 15)

	offset_range =(PresBBP >= 400.)
	BBP700_QC = BBP700_RM - BBP700_RM[offset_range].mean()

	BBP700_QC[BBP700_QC < 0.] = 0.

	return BBP700_QC

def CHL_QC(CHL):
	CHL_MF = ndimage.median_filter(CHL, size=5)
	CHL_QC = running_mean(CHL_MF, 7)
	CHL_QC[CHL_QC < 0.] = 0.
	return CHL_QC

def profile_shape(x, y):
	return x*y/np.max(y)


# Create functions for chl-specific IOP spectra

def PFT_calc(CHL, p1, p2, p3, p4):
	#p1, p2, p3 and p4 are relative contributions (0-1)
	# of various PFT to total absorption
	PFT_1 = p1*CHL     # 
	PFT_2 = p2*CHL
	PFT_3 = p3*CHL
	PFT_4 = p4*CHL
	return PFT_1, PFT_2, PFT_3, PFT_4

# Di Cicco et al., 2017 - both for PFT and PSC
def PFT_MED(CHL):
	x     = np.where(CHL>0, np.log10(CHL), 0.)
	#x = np.log10(CHL)
	MICRO = 0.0667*x**3 + 0.1939*x**2 + 0.2743*x + 0.2994
	NANO  =              -0.1740*x**2 - 0.0851*x + 0.4725
	PICO  = 1 - MICRO - NANO
	DIATOM = 0.0482*x**3 + 0.1877*x**2 + 0.2946*x + 0.2533
	DINOPH = MICRO - DIATOM
	CRYPT  = 0.0171*x**3 + 0.0667*x**2 + 0.1153*x + 0.0952
	GREEN  = (np.exp(-1.5780*x + 2.1841) + 22.6833 *x) ** (-1.)
	PROK   = 0.0664*x**3 + 0.1410*x**2 - 0.2097*x + 0.0979
	HAPT   = 1. - MICRO - CRYPT - GREEN - PROK
	
	PFT_1 = DIATOM * CHL
	PFT_2 = HAPT   * CHL
	PFT_3 = CRYPT  * CHL
	PFT_4 = GREEN  * CHL
	PFT_5 = PROK   * CHL
	PFT_6 = DINOPH  * CHL

	PSC_1 = MICRO  * CHL
	PSC_2 = NANO   * CHL
	PSC_3 = PICO   * CHL
	return PFT_1, PFT_2, PFT_3, PFT_4, PFT_5, PFT_6, PSC_1, PSC_2, PSC_3

	
def aNAP_Case1(CHL, Snap):

	aNAP = np.zeros((CHL.shape[0], wl.shape[0]))
	a440  = 0.0136*np.power(CHL, 0.615)

	for iwl in range(len(wl)):
		a_NAP = a440 * np.exp(-Snap * (wl[iwl] - 440.))
		aNAP[:,iwl] = a_NAP * CHL / np.max(CHL)

	return aNAP

def aCDOM_Case1(CHL, Scdom):

	aCDOM = np.zeros((CHL.shape[0], wl.shape[0]))
	a443   = 0.0316*np.power(CHL,0.63)

	for iwl in range(len(wl)):
		a_cdom = a443 * np.exp(-Scdom*(wl[iwl]-443.))
		aCDOM[:,iwl]  = a_cdom * CHL / np.max(CHL)   # tHIS IS GOING TO BE MODIFIED

	return aCDOM


def aCDOM_Case1_CDOM(CHL, CDOM, Scdom):

	aCDOM = np.zeros((CHL.shape[0], wl.shape[0]))
	a443   = 0.0316*np.power(CHL,0.63)

	CDOM_qc = CDOM_QC(CDOM)

	for iwl in range(len(wl)):
		a_cdom = a443 * np.exp(-Scdom*(wl[iwl]-443.))
		aCDOM[:,iwl]  = a_cdom * CDOM_qc / np.max(CDOM)   # tHIS IS GOING TO BE MODIFIED

	return aCDOM
