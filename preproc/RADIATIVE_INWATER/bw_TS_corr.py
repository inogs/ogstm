#!/bin/env python
# Xiaodong Zhang, Lianbo Hu, and Ming-Xia He (2009), Scatteirng by pure
# seawater: Effect of salinity, Optics Express, Vol. 17, No. 7, 5698-5710
#
# Main function of the module is betasw
from __future__ import print_function, division
import numpy as np

from configuration import *
from save_waterIOP import bw

Na  = 6.0221417930e23 # Avogadro's constant
Kbz = 1.3806503e-23   # Boltzmann constant
M0  = 18e-3           # Molecular weigth of water in kg/mol


def RInw(wv,Tc,S):
	'''
	Auxiliar function: Refractive index of seawater
	'''
	# Refractive index of air is from Ciddor (1996,Applied Optics)
	n_air = 1.0 + ( 5792105.0/(238.0185 - 1./(wv/1e3)*(wv/1e3) ) + 
		    167917.0/(57.362 - 1./(wv/1e3)*(wv/1e3)) )/1e8
	
	# Refractive index of seawater is from Quan and Fry (1994, Applied Optics)
	n0, n1, n2, n3, n4 = 1.31405, 1.779e-4, -1.05e-6, 1.6e-8, -2.02e-6
	n5, n6, n7, n8, n9 = 15.868, 0.01155, -0.00423, -4382, 1.1455e6

	nsw = n0 + (n1 + n2*Tc + n3*Tc*Tc)*S + n4*Tc*Tc + (n5 + n6*S + n7*Tc)/wv + n8/wv/wv + n9/wv/wv/wv # pure seawater
	nsw = nsw*n_air
	dnswds = (n1 + n2*Tc + n3*Tc*Tc + n6/wv)*n_air

	return nsw, dnswds

def BetaT(Tc,S):
	'''
	Auxiliar function: isothermal compressibility
	'''
	# pure water secant bulk Millero (1980, Deep-sea Research)
	kw = 19652.21 + 148.4206*Tc - 2.327105*Tc*Tc + 1.360477e-2*Tc*Tc*Tc - 5.155288e-5*Tc*Tc*Tc*Tc;
	Btw_cal = 1./kw;

	# isothermal compressibility from Kell sound measurement in pure water
	# Btw = (50.88630 + 0.717582*Tc + 0.7819867e-3*Tc*Tc + 31.62214e-6*Tc*Tc*Tc - 0.1323594e-6*Tc*Tc*Tc*Tc + 0.634575e-9*Tc*Tc*Tc*Tc*Tc)/(1 + 21.65928e-3*Tc)*1e-6

	# seawater secant bulk
	a0 = 54.6746  - 0.603459*Tc  + 1.09987e-2*Tc*Tc - 6.167e-5*Tc*Tc*Tc
	b0 = 7.944e-2 + 1.6483e-2*Tc - 5.3009e-4*Tc*Tc

	Ks = kw + a0*S + b0*S**1.5

	# calculate seawater isothermal compressibility from the secant bulk
	return 1./Ks*1e-5 # unit is pa

def rhou_sw(Tc,S):
	'''
	Auxiliar function: seawater density in Kg/m^3, from UNESCO,38,1981
	'''
	a0 = 8.24493e-1;   a1 = -4.0899e-3;  a2 = 7.6438e-5;   a3 = -8.2467e-7; a4 = 5.3875e-9;
	a5 = -5.72466e-3;  a6 = 1.0227e-4;   a7 = -1.6546e-6;  a8 = 4.8314e-4;
	b0 = 999.842594;   b1 = 6.793952e-2; b2 = -9.09529e-3; b3 = 1.001685e-4;
	b4 = -1.120083e-6; b5 = 6.536332e-9;

	# density for pure water 
	density_w = b0 + b1*Tc + b2*Tc*Tc + b3*Tc*Tc*Tc + b4*Tc*Tc*Tc*Tc + b5*Tc*Tc*Tc*Tc*Tc
	# density for pure seawater
	density_sw = density_w +( 
			(a0 + a1*Tc + a2*Tc*Tc + a3*Tc*Tc*Tc + a4*Tc*Tc*Tc*Tc)*S +
			(a5 + a6*Tc + a7*Tc*Tc)*S**1.5 + a8*S*S
		)

	return density_sw

def dlnasw_ds(Tc,S):
	'''
	Auxiliar function: 
	water activity data of seawater is from Millero and Leung (1976,American
	Journal of Science,276,1035-1077). Table 19 was reproduced using
	Eqs.(14,22,23,88,107) then were fitted to polynominal equation.
	dlnawds is partial derivative of natural logarithm of water activity
	w.r.t.salinity
	'''
#	lnaw = (-1.64555e-6 - 1.34779e-7*Tc + 1.85392e-9*Tc*Tc - 1.40702e-11*Tc*Tc*Tc) +
#		   (-5.58651e-4 + 2.40452e-7*Tc - 3.12165e-9*Tc*Tc + 2.40808e-11*Tc*Tc*Tc)*S +
#		   ( 1.79613e-5 - 9.9422e-8*Tc  + 2.08919e-9*Tc*Tc - 1.39872e-11*Tc*Tc*Tc)*S**1.5 +
#		   (-2.31065e-6 - 1.37674e-9*Tc - 1.93316e-11*Tc*Tc)*S*S
	return (-5.58651e-4 + 2.40452e-7*Tc - 3.12165e-9*Tc*Tc + 2.40808e-11*Tc*Tc*Tc) + 1.5*( 1.79613e-5 - 9.9422e-8*Tc  + 2.08919e-9*Tc*Tc - 1.39872e-11*Tc*Tc*Tc)*np.sqrt(S) + 2*(-2.31065e-6 - 1.37674e-9*Tc - 1.93316e-11*Tc*Tc)*S

def PMH(n_wat):
	'''
	Auxiliar function: density derivative of refractive index from PMH model
	'''
	n_wat2 = n_wat*n_wat
	return (n_wat2 - 1.)*( 1. + 2./3.*(n_wat2 + 2)*(n_wat/3. - 1./3./n_wat)**2 )

def betasw(wv,Tc,S,theta,delta=0.039):
	'''
	Xiaodong Zhang, Lianbo Hu, and Ming-Xia He (2009), Scatteirng by pure
	seawater: Effect of salinity, Optics Express, Vol. 17, No. 7, 5698-5710

	Inputs:
		> wv (nm):  wavelength
		> Tc:       temperature in degree Celsius, must be a scalar
		> S:        salinity, must be scalar
		> theta:    angles for volume scattering
		> delta:    depolarization ratio, if not provided, default = 0.039 will be used.

	Outputs:
		> bsw:      total scattering coefficient.
		> beta90sw: volume scattering at 90 degree.
		> betasw:   volume scattering at angles defined by theta.

	For backscattering coefficients, divide total scattering by 2

	Xiaodong Zhang, March 10, 2009 (Adapted to python)
	'''
	Tk  = Tc + 273.15       # Absolute tempearture
	rad = np.deg2rad(theta) # Angle in radian as a colum variable

	# nsw:  absolute refractive index of seawater
	# dnds: partial derivative of seawater refractive index w.r.t. salinity
	nsw, dnds = RInw(wv,Tc,S)

	# Isothermal compressibility is from 
	# Lepple & Millero (1971,Deep Sea-Research), pages 10-11
	# The error ~ +/-0.004e-6 bar^-1
	IsoComp = BetaT(Tc,S)

	# density of water and seawater,unit is Kg/m^3, from UNESCO,38,1981
	density_sw = rhou_sw(Tc,S)

	# water activity data of seawater is from Millero and Leung (1976,American
	# Journal of Science,276,1035-1077). Table 19 was reproduced using
	# Eq.(14,22,23,88,107) then were fitted to polynominal equation.
	# dlnawds is partial derivative of natural logarithm of water activity
	# w.r.t.salinity
	dlnawds = dlnasw_ds(Tc,S)

	# density derivative of refractive index from PMH model
	DFRI = PMH(nsw) # PMH model

	# volume scattering at 90 degree due to the density fluctuation
	beta_df = np.pi*np.pi/2.*( (wv*1e-9)**(-4.) )*Kbz*Tk*IsoComp*DFRI*DFRI*(6. + 6.*delta)/(6. - 7.*delta)

	# volume scattering at 90 degree due to the concentration fluctuation
	flu_con = S*M0*dnds*dnds/density_sw/(-dlnawds)/Na
	beta_cf = 2.*np.pi*np.pi*( (wv*1e-9)**(-4) )*nsw*nsw*(flu_con)*(6. + 6.*delta)/(6. - 7.*delta)

	# total volume scattering at 90 degree
	beta90sw = beta_df+beta_cf
	bsw = 8.*np.pi/3.*beta90sw*(2. + delta)/(1 + delta)
	betasw = np.zeros((len(theta),len(wv)))
	for iwv in range(len(wv)):
		betasw[:,iwv] = beta90sw[iwv]*( 1. + np.cos(rad)*np.cos(rad)*(1. - delta)/(1. + delta) )

	return bsw, beta90sw, betasw


def bsw(wv,Tc,S,delta=0.039):
	'''
	Xiaodong Zhang, Lianbo Hu, and Ming-Xia He (2009), Scatteirng by pure
	seawater: Effect of salinity, Optics Express, Vol. 17, No. 7, 5698-5710

	Inputs:
		> wv (nm):  wavelength
		> Tc:       temperature in degree Celsius, must be a scalar
		> S:        salinity, must be scalar

		> delta:    depolarization ratio, if not provided, default = 0.039 will be used.

	Outputs:
		> bsw:      total scattering coefficient.

	For backscattering coefficients, divide total scattering by 2

	Xiaodong Zhang, March 10, 2009 (Adapted to python)
	'''
	Tk  = Tc + 273.15       # Absolute tempearture

	# nsw:  absolute refractive index of seawater
	# dnds: partial derivative of seawater refractive index w.r.t. salinity
	nsw, dnds = RInw(wv,Tc,S)

	# Isothermal compressibility is from 
	# Lepple & Millero (1971,Deep Sea-Research), pages 10-11
	# The error ~ +/-0.004e-6 bar^-1
	IsoComp = BetaT(Tc,S)

	# density of water and seawater,unit is Kg/m^3, from UNESCO,38,1981
	density_sw = rhou_sw(Tc,S)

	# water activity data of seawater is from Millero and Leung (1976,American
	# Journal of Science,276,1035-1077). Table 19 was reproduced using
	# Eq.(14,22,23,88,107) then were fitted to polynominal equation.
	# dlnawds is partial derivative of natural logarithm of water activity
	# w.r.t.salinity
	dlnawds = dlnasw_ds(Tc,S)

	# density derivative of refractive index from PMH model
	DFRI = PMH(nsw) # PMH model

	# volume scattering at 90 degree due to the density fluctuation
	beta_df = np.pi*np.pi/2.*( (wv*1e-9)**(-4.) )*Kbz*Tk*IsoComp*DFRI*DFRI*(6. + 6.*delta)/(6. - 7.*delta)

	# volume scattering at 90 degree due to the concentration fluctuation
	flu_con = S*M0*dnds*dnds/density_sw/(-dlnawds)/Na
	beta_cf = 2.*np.pi*np.pi*( (wv*1e-9)**(-4) )*nsw*nsw*(flu_con)*(6. + 6.*delta)/(6. - 7.*delta)

	# total volume scattering at 90 degree
	beta90sw = beta_df+beta_cf
	bsw = 8.*np.pi/3.*beta90sw*(2. + delta)/(1 + delta)
	
	return bsw


def bw_TS_corr(T, S):

    bw_TS = np.zeros((wl.shape[0], T.shape[0]))

    for depth in range(len(T)):

        bw_TS[:,depth] = bsw(wl, T[depth] , S[depth], delta=0.039)

    return bw_TS


def bw_NO_corr(T, S):

    bw_NO = np.zeros((wl.shape[0], T.shape[0]))

    for depth in range(len(T)):

        bw_NO[:,depth] = bw #bsw(wl, T[depth] , S[depth], delta=0.039)

    return bw_NO


