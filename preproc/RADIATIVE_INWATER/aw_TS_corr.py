#!/bin/env python

from __future__ import print_function
import numpy as np
from scipy import interpolate

from configuration import *
from save_waterIOP import aw_OASIM, aw_LIT, aw_MASON

'''
With this function you correct for the T-S contribution to
pure water absorption .
For this we use Sullivan et al. 2006 (AO) 
'''

mydtype= np.dtype([
          ('wl'         ,      np.int),
          ('phi_t'      ,    np.float), # not np.float32, it performs a minor change
          ('sgm_phi_t'  ,    np.float),
          ('phi_s_c'    ,    np.float),
          ('sgm_phi_s_c',    np.float),
          ('phi_s_a'    ,    np.float),
          ('sgm_phi_s_a',    np.float)] )

INPUT = np.loadtxt('Sullivan_T_chart.txt',dtype=mydtype, delimiter=" ",skiprows=1, ndmin=1)

wavelengths = INPUT['wl']
phi_t       = INPUT['phi_t']   # Temperature parameter
phi_s_a     = INPUT['phi_s_a'] # Salinity parameter for absorption. c stands for beam attenuation, we don't need that here

# Now we interpolate Sullivan's phi_T values to our model wavelengths - imported from configuration.py

a_phi_t = interpolate.interp1d(wavelengths, phi_t,   kind='linear', fill_value='extrapolate')(wl)
a_phi_s = interpolate.interp1d(wavelengths, phi_s_a, kind='linear', fill_value='extrapolate')(wl)

# We are going to test 2 different methods to estimate aw: Pope and Fry vs. Mason's.
# Each has a different Tref, where measurements were taken

Tref_OASIM = 22.     # OASIM:  Tref = 22 deg C+/- 1 
Tref_LIT   = 22.     # Lite.:  Tref = 22 deg C+/- 1 
Tref_MASON = 23.     # Mason:  Tref = 23 deg C +/-  0.5 deg C         

def aw_TS_corr(T, S, model='MASON'):

    aw_TS = np.zeros((a_phi_t.shape[0], T.shape[0]))
    
    dT = Tref_MASON - T if model=='MASON' else Tref_OASIM - T
    dS = S

    aw = aw_MASON if model=='MASON' else aw_OASIM 

    for depth in range(len(T)):

        aw_TS[:,depth] = aw - a_phi_t*dT[depth] + a_phi_s*dS[depth]

    return aw_TS

def aw_NO_corr(T, S, model='MASON'):

    aw_NO = np.zeros((a_phi_t.shape[0], T.shape[0]))
    
    dT = Tref_MASON - T if model=='MASON' else Tref_OASIM - T
    dS = S

    aw = aw_MASON if model=='MASON' else aw_OASIM 

    for depth in range(len(T)):

        aw_NO[:,depth] = aw  #- a_phi_t*dT[depth] + a_phi_s*dS[depth]

    return aw_NO
