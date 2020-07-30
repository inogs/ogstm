#!/bin/env python

from __future__ import print_function
import numpy as np
from scipy import interpolate

from configuration import *

'''With this function you correct for the T-S contribution to
pure water absorption .
For this we use Sullivan et al. 2006 (AO) '''

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

a_phi_t = interpolate.interp1d(wl_phi, phi_t,   kind='linear', fill_value='extrapolate')(wl)
a_phi_s = interpolate.interp1d(wl_phi, phi_s_a, kind='linear', fill_value='extrapolate')(wl)

# Initialize an array of ones      
init = np.ones(len(wl))

# Final array

# We are going to test 2 different methods to estimate aw: Pope and Fry vs. Mason's.
# Each has a different Tref, where measurements were taken

# OASIM:  Tref = 22 °C +/- 1 °C
# Mason:  Tref = 23 °C +/-  0.5°C

Tref_OASIM = 22.    ; dT_OASIM = T - Tref_OASIM   ; dS = S
Tref_MASON = 23.    ; dT_MASON = T - Tref_MASON   ; dS = S

a_TS_OASIM = a - init*[i*dT_OASIM for i in a_phi_t] + init*[j*dS for j in a_phi_s]
