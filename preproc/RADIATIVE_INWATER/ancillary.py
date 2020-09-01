#!/bin/env python

import netCDF4 as NC4
import numpy as np
import os
from scipy import ndimage
from scipy.optimize import curve_fit

from aPHY_read import aPFT1, aPFT2, aPFT3, aPFT4, aPFT5, aPFT6
from aw_TS_corr import *
from bw_TS_corr import rhou_sw
from configuration import wl
from Q_read import Qo_morel, SQn_morel



''' Here you put all the functions, also for the IOPs , T-S corrections, etc. '''

def str2bool(v):
  return True if v.lower() in ("yes", "true", "t", "1") else False

def INT_wl(xa, xb, ya, yb, xINT):
    yINT = ya*(xb-xINT)/(xb-xa) + yb*(xa-xINT)/(xa-xb)
    return yINT


def write_abw25(wl, aw, bw, fname='bcs/abw25.dat'):

	file  = open(fname, 'w')
	for idepth in range(aw.shape[1]):
		for iwl in range(aw.shape[0]):
			file.write('%5d%15.4f%10.4f\n'%(wl[iwl], aw[iwl, idepth], bw[iwl, idepth]))
	file.close()

	return 


def write_acbc25(wl, aPFT_TOT, bp_TOT, bbp_TOT, fname='bcs/acbc25.dat'):

	file  = open(fname, 'w')
	for idepth in range(aPFT_TOT.shape[0]):
		for iwl in range(aPFT_TOT.shape[1]):
			file.write('%5d%15.8f%15.8f%15.8f\n'%(wl[iwl], aPFT_TOT[idepth, iwl], bp_TOT[idepth, iwl], bbp_TOT[idepth, iwl]))
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

def save_matchup(ncfile, Pres380, Pres412, Pres490, Ed380_float, Ed412_float, Ed490_float, Ed380_model, Ed412_model, Ed490_model, timestr):

	#if not os.path.isdir('MATCHUP') : os.mkdir('MATCHUP')

	modelfile = 'MATCHUP/' + ncfile
	ncmodel   = NC4.Dataset(modelfile,"w");
				
	ncdepth380 = ncmodel.createDimension('depth380',     len(Pres380));
	ncdepth412 = ncmodel.createDimension('depth412',     len(Pres412));
	ncdepth490 = ncmodel.createDimension('depth490',     len(Pres490));

	ncwave  = ncmodel.createDimension('type', 2);

	setattr(ncmodel, 'time', timestr);
	
	ncDepth380 = ncmodel.createVariable('depth380', 'f', ('depth380')); 
	setattr(ncDepth380, 'unit',  '[m]' );
	ncDepth380[:] = Pres380

	ncDepth412 = ncmodel.createVariable('depth412', 'f', ('depth412')); 
	setattr(ncDepth412, 'unit',  '[m]' );
	ncDepth412[:] = Pres412

	ncDepth490 = ncmodel.createVariable('depth490', 'f', ('depth490')); 
	setattr(ncDepth490, 'unit',  '[m]' );
	ncDepth490[:] = Pres490
	
	ncEd380 = ncmodel.createVariable('Ed_380', 'f', ('depth380', 'type'));
	setattr(ncEd380, 'missing_value',-1.0 );     
	setattr(ncEd380, 'long_name',  'Downward irradiance at 380 nm ' );     
	setattr(ncEd380, 'unit',  '[W/m2/nm]' );
	ncEd380[:] = np.vstack((Ed380_float, Ed380_model)).T

	ncEd412 = ncmodel.createVariable('Ed_412', 'f', ('depth412', 'type'));
	setattr(ncEd412, 'missing_value',-1.0 );     
	setattr(ncEd412, 'long_name',  'Downward irradiance at 412 nm ' );     
	setattr(ncEd412, 'unit',  '[W/m2/nm]' );
	ncEd412[:] = np.vstack((Ed412_float, Ed412_model)).T

	ncEd490 = ncmodel.createVariable('Ed_490', 'f', ('depth490', 'type'));
	setattr(ncEd490, 'missing_value',-1.0 );     
	setattr(ncEd490, 'long_name',  'Downward irradiance at 490 nm ' );     
	setattr(ncEd490, 'unit',  '[W/m2/nm]' );
	ncEd490[:] = np.vstack((Ed490_float, Ed490_model)).T
	
	ncmodel.close()

	return ncmodel


def save_Kd(ncfile, Kd_model, Kd_float, wl_Kd, timestr):

	#if not os.path.isdir('KD') : os.mkdir('KD')

	modelfile = 'KD/' + ncfile
	ncmodel   = NC4.Dataset(modelfile,"w");
				
	ncdepth = ncmodel.createDimension('depth',     1);
	ncwave  = ncmodel.createDimension('wavelength', len(wl_Kd));

	setattr(ncmodel, 'time', timestr);
	
	ncWL = ncmodel.createVariable('wavelength', 'f', ('wavelength')); 
	setattr(ncWL, 'unit',  '[nm]' );
	ncWL[:] = wl_Kd
	
	ncKdM = ncmodel.createVariable('Kd_model', 'f', ('depth', 'wavelength'));
	setattr(ncKdM, 'missing_value',-1.0 );     
	setattr(ncKdM, 'long_name',  'Diffuse attenuation coefficient for downward planar irradiance' );     
	setattr(ncKdM, 'unit',  '[1/m]' );

	ncKdM[:] = Kd_model

	ncKdF = ncmodel.createVariable('Kd_float', 'f', ('depth', 'wavelength'));
	setattr(ncKdF, 'missing_value',-1.0 );     
	setattr(ncKdF, 'long_name',  'Diffuse attenuation coefficient for downward planar irradiance' );     
	setattr(ncKdF, 'unit',  '[1/m]' );

	ncKdF[:] = Kd_float
	
	ncmodel.close()

	return ncmodel

def save_reflectance(ncfile, wl_RRS, Rrs, timestr):

	#if not os.path.isdir('RRS') : os.mkdir('RRS')

	modelfile = 'RRS/' + ncfile
	ncmodel   = NC4.Dataset(modelfile,"w");
				
	ncdepth = ncmodel.createDimension('depth',     1);
	ncwave  = ncmodel.createDimension('wavelength', len(wl_RRS));

	setattr(ncmodel, 'time', timestr);
	
	ncWL = ncmodel.createVariable('wavelength', 'f', ('wavelength')); 
	setattr(ncWL, 'unit',  '[nm]' );
	ncWL[:] = wl_RRS
	
	ncRrs = ncmodel.createVariable('Rrs', 'f', ('depth', 'wavelength'));
	setattr(ncRrs, 'missing_value',-1.0 );     
	setattr(ncRrs, 'long_name',  'Remote sensing reflectance ' );     
	setattr(ncRrs, 'unit',  '[1/sr]' );

	ncRrs[:] = Rrs
	
	ncmodel.close()

	return ncmodel

# Calculate Kw based on aw and bw data from our models

def calc_Kw380(aw380, bw380):
	Kw_380 = aw380 + 0.5 * bw380
	return Kw_380

# Calculating euphotic and MLD depth ranges for Kd calculation

def euphotic(Pres, Ed):
	EUPH               = [x for x in range(len(Pres)) if Ed[x]/Ed[0] > 0.37]
	ind                = EUPH[-1]
	zeu                = Pres[ind]
	Pres_Euph          = Pres[EUPH]
	Ed_Euph            = Ed[EUPH]
	return zeu, Pres_Euph, Ed_Euph 

def MLD_calc(Pres, T, S, Ed):

	density  = rhou_sw(T, S)           # Potential density

	ind      = np.argmin(np.abs(Pres - 10.))  # Find the index where depth is closest to 10 m.
	MLD_ind  = np.argmin(np.abs(density[ind] - density) - 0.03)  # Find the index where the potential density is closest to 0.03

	zMLD     = Pres[MLD_ind]           # MLD depth
	Pres_MLD = Pres[0:MLD_ind]         # MLD depth range
	Ed_MLD   = Ed[  0:MLD_ind]           # Ed at the MLD depth range

	return zMLD, Pres_MLD, Ed_MLD

def calc_Kd(Pres, Ed):

	success = True

	# Optical depth  range (DR) calculation
	zmax, PresDR, Ed_DR = euphotic(Pres, Ed) 

	if len(PresDR) < 5 or PresDR[1] - PresDR[0] > 9. or PresDR[4] > 15.:
		success = False
		Kd      = np.nan

	func = lambda z, Kd: np.exp(-Kd*z)

	if success:
		popt, pcov = curve_fit(func, PresDR-PresDR[0], Ed_DR/Ed_DR[0], p0=0.1)
		Kd         = popt[0]

	return Kd, success  

def calc_Kbio_380(Ed_380, Pres, T, S, depth_type, aw380, bw380, Kw_type):

	success = True

	# Euphotic or MLD depth range (DR) calculation
	zmax, PresDR, Ed_DR = euphotic(Pres, Ed_380) if depth_type == 'EUPH' else MLD_calc(Pres, T, S, Ed_380)

	if len(PresDR) < 5 or PresDR[1] - PresDR[0] > 9. or PresDR[4] > 15.:
		success = False
		Kd_380  = np.nan

	func = lambda z, Kd: np.exp(-Kd*z)

	if success:
		popt, pcov = curve_fit(func, PresDR-PresDR[0], Ed_DR/Ed_DR[0], p0=0.1)
		Kd_380 = popt[0]

	Kw_380 = np.ones((Pres.shape[0]))*0.01510 if Kw_type == 'MOREL' else calc_Kw380(aw380, bw380)  # According to Morel and Maritorena (2001)

	Kbio_380 = Kd_380 - Kw_380

	if np.any(Kbio_380 < 0):
		success = False

	return Kbio_380, success   # needs to be a returned array
 
# BGC-Argo profile QC procedures

def CDOM_QC(CDOM):
	CDOM_MF = ndimage.median_filter(CDOM, size=5)
	CDOM_QC = running_mean(CDOM_MF, 7)
	CDOM_QC[CDOM_QC < 0.] = 0.
	
	return CDOM_QC

def BBP700_QC(PresBBP, BBP700):
	BBP700_MF = ndimage.median_filter(BBP700, size=11)
	BBP700_QC = running_mean(BBP700_MF, 15)

	BBP700_QC[BBP700_QC < 0.] = 0.

	return BBP700_QC

def CHL_QC(CHL):  # For now we are not using it
	CHL_MF = ndimage.median_filter(CHL, size=5)
	CHL_QC = running_mean(CHL_MF, 7)
	CHL_QC[CHL_QC < 0.] = 0.
	return CHL_QC

def profile_shape(x, y):
	return x*y/np.max(y)


# Create functions for chl-specific IOP spectra

# Di Cicco et al., 2017 - both for PFT and PSC
def PFT_MED(CHL):
# This algorithm works only if CHL is in range between 0.02 and 5.52 mg/m3
# And it was validated only for the first 50 meters
	x      = np.where(CHL>0.02, np.log10(CHL), np.log10(0.02))
	x[CHL>5.52] = np.log10(5.52)
	#x = np.log10(CHL)
	MICRO  = 0.0667*x**3 + 0.1939*x**2 + 0.2743*x + 0.2994
	NANO   =              -0.1740*x**2 - 0.0851*x + 0.4725
	PICO   = 1 - MICRO - NANO
	DIATOM = 0.0482*x**3 + 0.1877*x**2 + 0.2946*x + 0.2533
	DINOPH = MICRO - DIATOM
	CRYPT  = 0.0171*x**3 + 0.0667*x**2 + 0.1153*x + 0.0952
	GREEN  = (np.exp(-1.5780*x + 2.1841) + 22.6833 *x) ** (-1.)
	PROK   = 0.0664*x**3 + 0.1410*x**2 - 0.2097*x + 0.0979
	HAPT   = 1. - MICRO - CRYPT - GREEN - PROK

	PFT_1  = DIATOM  * CHL
	PFT_2  = DINOPH  * CHL
	PFT_3  = CRYPT   * CHL
	PFT_4  = GREEN   * CHL
	PFT_5  = PROK    * CHL
	PFT_6  = HAPT    * CHL
 
	PSC_1  = MICRO   * CHL
	PSC_2  = NANO    * CHL
	PSC_3  = PICO    * CHL
 
	aPFT_TOT = np.zeros((CHL.shape[0], wl.shape[0]))

	for idepth in range(CHL.shape[0]):

		aPFT_TOT[idepth, :] = aPFT1*PFT_1[idepth] + aPFT2*PFT_2[idepth] + aPFT3*PFT_3[idepth] + aPFT4*PFT_4[idepth] + aPFT5*PFT_5[idepth] + aPFT6*PFT_6[idepth]

	return aPFT_TOT #PFT_1, PFT_2, PFT_3, PFT_4, PFT_5, PFT_6, PSC_1, PSC_2, PSC_3


def aNAP_Babin(X, a443, Snap):

	# X is CHL or BBP700

	aNAP = np.zeros((X.shape[0], wl.shape[0]))

	for iwl in range(len(wl)):
		a_NAP = a443 * np.exp(-Snap * (wl[iwl] - 443.))
		aNAP[:,iwl] = a_NAP * X / np.max(X)

	return aNAP


def aNAP_Case1(CHL, X, Snap):

	aNAP = np.zeros((CHL.shape[0], wl.shape[0]))
	a440  = 0.0136*np.power(CHL, 0.615)

	for iwl in range(len(wl)):
		a_NAP = a440 * np.exp(-Snap * (wl[iwl] - 440.))
		aNAP[:,iwl] = a_NAP * X / np.max(X)

	return aNAP

def aCDOM_Case1(CHL, Scdom):

	aCDOM = np.zeros((CHL.shape[0], wl.shape[0]))
	a443   = 0.0316*np.power(CHL,0.63)

	for iwl in range(len(wl)):
		a_cdom = a443 * np.exp(-Scdom*(wl[iwl]-443.))
		aCDOM[:,iwl]  = a_cdom * CHL / np.max(CHL)   

	return aCDOM


def aCDOM_Case1_CDOM(CHL, CDOM_qc, Scdom):

	aCDOM  = np.zeros((CHL.shape[0], wl.shape[0]))
	a443   = 0.0316*np.power(CHL,0.63)

	for iwl in range(len(wl)):
		a_cdom = a443 * np.exp(-Scdom*(wl[iwl]-443.))
		aCDOM[:,iwl]  = a_cdom * CDOM_qc / np.max(CDOM_qc) 

	return aCDOM

def aCDOM_Kbio(CDOM_int, Scdom, Kbio_380):

	aCDOM  = np.zeros((CDOM_int.shape[0], wl.shape[0]))      # CDOM interpolated to PresCHL
	a380   = Kbio_380#*np.ones((CDOM_int.shape[0]))

	for iwl in range(len(wl)):
		a_cdom = a380 * np.exp(-Scdom*(wl[iwl]-380.))
		aCDOM[:,iwl]  = a_cdom * CDOM_int / np.max(CDOM_int) 

	return aCDOM

def aDG_CORR(aDG, aw):  # detritus (NAP) and gelbstoff (CDOM) correction of the spectra

	aDG_CORR = np.zeros((aDG.shape[0], wl.shape[0]))  

	T = np.zeros((aDG.shape[0],))
	S = np.zeros((aDG.shape[0],))

	aw_LIT = aw_NO_corr(T, S, model='LIT')

	for iwl in range(len(wl)):
		aDG_CORR[:,iwl] = aDG[:,iwl] + aw_LIT[iwl,:] - aw[iwl,:]
		
	return aDG_CORR

def bp_Case1(CHL, ratio):
	bp = np.zeros((CHL.shape[0], wl.shape[0]))

	for idepth in range(CHL.shape[0]):
		# Assigning the nu function for CHL < 0.02 to calculate it as if CHL were 0.02 mg m-3
		# In Morel et al., 2002 there is only the range between 0.02 and 2 mg m-3 and for above 2 mg m-3
		nu = 0.5 * (np.log10(CHL[idepth]) - 0.3) if CHL[idepth] > 0.02 else (np.log10(0.02) - 0.3)
		if CHL[idepth] > 2. : nu = 0.

		bp[idepth,:] = 0.416*np.power(CHL[idepth], 0.766)*np.power((wl/550.), nu)  * CHL[idepth] / np.max(CHL)  

	bbp = bp * ratio

	return bp, bbp

def bp_Case1_bbp(CHL, bbp700_int, ratio):
	bp = np.zeros((CHL.shape[0], wl.shape[0]))

	for idepth in range(CHL.shape[0]):
		# Assigning the nu function for CHL < 0.02 to calculate it as if CHL were 0.02 mg m-3
		# In Morel et al., 2002 there is only the range between 0.02 and 2 mg m-3 and for above 2 mg m-3
		nu = 0.5 * (np.log10(CHL[idepth]) - 0.3) if CHL[idepth] > 0.02 else (np.log10(0.02) - 0.3)
		if CHL[idepth] > 2. : nu = 0.

		bp[idepth,:] = 0.416*np.power(CHL[idepth], 0.766)*np.power((wl/550.), nu)  * bbp700_int[idepth] / np.max(bbp700_int)  

	bbp = bp * ratio

	return bp, bbp

def bbp_frombbp700(bbp700_int, slope, ratio):
	bbp = np.zeros((bbp700_int.shape[0], wl.shape[0]))

	for iwl in range(wl.shape[0]):
		bbp[:, iwl] = bbp700_int * np.power((wl[iwl]/700.), -slope) * bbp700_int/ np.max(bbp700_int)  

	bp = bbp / ratio

	return bp, bbp


# Compute Q

def Q_aas(solz):
	Q = 5.33 * np.exp(-0.45 * np.cos(np.deg2rad(solz)))
	return Q



def Q_morel(solz, PresCHL, CHL, wl, depth_type):

	if depth_type == 'SURF':
		Qo  = Qo_morel( wl, CHL[0])[0]
		Sqn = SQn_morel(wl, CHL[0])[0]
	if depth_type == '10m':
		Qo  = Qo_morel( wl, np.mean(CHL[PresCHL<10.]))[0]
		Sqn = SQn_morel(wl, np.mean(CHL[PresCHL<10.]))[0]

	Q = Qo + Sqn * (1 - np.cos(np.deg2rad(solz)))

	return Q




