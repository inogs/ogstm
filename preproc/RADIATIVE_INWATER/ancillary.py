#!/bin/env python

import netCDF4 as NC4
import numpy as np

''' Here you put all the functions, also for the IOPs , T-S corrections, etc. '''

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

# Create functions for chl-specific IOP spectra

# Wavelengths of the Radiative Transfer Model
lam = np.array([  250.0, 325.0, 350.0, 375.0, 400.0, 425.0,   450.0,  475.0, 500.0 , 525.0, 550.0, 575.0, 600.0, 625.0, 650.0, 675.0, 700.0, 725.0, 775.0, 
				  850.0, 950.0, 1050.0, 1150.0, 1250.0, 1350.0, 1450.0, 1550.0,  1650.0, 1750.0, 1900.0, 2200.0, 2900.0, 3700.0 , 4000.0]) 
lam = lam.reshape(lam.shape[0],1)

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
	

def NAP_abs(CHL, Snap, a443):
	a_NAP = a443 * np.exp(-Snap * (lam - 443.))
	a_NAP = a_NAP.reshape((1, a_NAP.shape[0]))
	CHL = CHL.reshape(CHL.shape[0], 1)
	aNAP = a_NAP * CHL / np.max(CHL)
	return aNAP

def CDOM_abs(CHL, Scdom, a440):
	a_cdom = a440 * np.exp(-Scdom*(lam-440.))
	a_cdom = a_cdom.reshape((1, a_cdom.shape[0]))
	CHL = CHL.reshape(CHL.shape[0], 1)
	aCDOM = a_cdom * CHL / np.max(CHL)
	return aCDOM


def profile_shape(x, y):
	
	return x*y/np.max(y)
