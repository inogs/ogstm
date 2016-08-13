#! /usr/bin/python

#AUTHOR PL 15.X.2013

# LOAD PACKAGES

import os,sys

import collections

import numpy as np

from mydtype import *

import scipy.io.netcdf as NC

import pickle

import do_big_ave as dba
import do_big_ave_phys as dbap
import do_big_ave_dia as dbad
import plot_hovmoeller_PTSW as plot_PTSW
import plot_hovmoeller_LDNCCC as plot_LDNCCC
import plot_hovmoeller_LDNC   as plot_LDNC  
import plot_hovmoeller_CCCCBBBB   as plot_CCCCBBBB  
import plot_RL_NP as plot_TIL
import plot_RL_NP_movie as plot_TIL_mov


# MAIN PROGRAM

TEST_LIST=np.loadtxt('TEST_LIST.dat', dtype=test_conf,skiprows=1,ndmin=1)

for test in TEST_LIST:

    print test['Dir']

#   dba.do_big_ave(test)
#   dbap.do_big_ave_phys(test)
#   dbad.do_big_ave_dia(test)
#   plot_PTSW.plot_hovmoeller_PTSW(test)
#   plot_LDNCCC.plot_hovmoeller_LDNCCC(test)
#   plot_LDNC.plot_hovmoeller_LDNC(test)
#   plot_CCCCBBBB.plot_hovmoeller_CCCCBBBB(test)
    plot_TIL_mov.plot_RL_NP_movie(test)

