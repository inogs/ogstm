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
import plot_hovmoeller_PTSW as plot_PTSW


# MAIN PROGRAM

TEST_LIST=np.loadtxt('TEST_LIST.dat', dtype=test_conf,skiprows=1,ndmin=1)

for test in TEST_LIST:

    print test['Dir']

#   dba.do_big_ave(test)
    dbap.do_big_ave_phys(test)
    plot_PTSW.plot_hovmoeller_PTSW(test)

