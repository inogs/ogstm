#! /usr/bin/python

#AUTHOR PL 15.X.2013

# LOAD PACKAGES

import os,sys

import collections

import numpy as np

from mydtype import *

import scipy.io.netcdf as NC

import pickle

import create_meshmask_nc as c_mask 

import create_Dom_Dec as c_dom 

import create_extinction_nc as c_ext

import create_init_nc as c_init

import create_forcings_nc as c_for

import create_bc_nc as c_bc

import deploy_code as d_code

import create_events as c_events

import create_da_nc as DA

# MAIN PROGRAM

TEST_LIST=np.loadtxt('TEST_LIST.dat', dtype=test_conf,skiprows=1,ndmin=1)

for test in TEST_LIST:

    print test['Dir']

    DA.create_dataset(test)

    c_dom.create_Dom_Dec(test)
    
    c_mask.create_meshmask_nc(test)
    
    c_for.create_forcings_nc(test)
    
    c_ext.create_extinction_nc(test)
    
    c_bc.create_bc_nc(test)
    
    c_init.create_init_nc(test)

    d_code.deploy_code(test)

    c_events.create_events(test)
