#! /usr/bin/env python

#AUTHOR PL 15.X.2013

# LOAD PACKAGES

import os,sys

import collections

import numpy as np

from mydtype import *

import scipy.io.netcdf as NC

import pickle

import do_big_ave as dba
import do_big_ave_HF as dbaHF
import do_big_ave_phys as dbap
#import do_big_ave_dia as dbad
#import do_big_ave_flux as dbaf
#import plot_hovmoeller_PTSW as plot_PTSW
#import plot_hovmoeller_LDNCCC as plot_LDNCCC
#import plot_hovmoeller_LDNC   as plot_LDNC  
#import plot_hovmoeller_CCCCBBBB   as plot_CCCCBBBB  
import plot_hovmoeller_C_CF as plot_C_CF
#import plot_COMPARE_DCM as plot_CMP_DCM
#import plot_hovmoeller_BUDGET as plt_BUD

try:
    from mpi4py import MPI
    comm  = MPI.COMM_WORLD
    rank  = comm.Get_rank()
    nranks =comm.size
    isParallel = True
except:
    rank   = 0
    nranks = 1
    isParallel = False

# MAIN PROGRAM

#LIST_POSTPROC = ['POSTPROC_REF','POSTPROC_Kd_1e-04', 'POSTPROC_THETA_DIA_0_5', 'POSTPROC_p_srsX2', 'POSTPROC_alpha-10p', 'POSTPROC_alpha-20p']
LIST_POSTPROC  = ['POSTPROC']
POSTPROC_DIR   = 'POSTPROC'

TEST_LIST=np.loadtxt('TEST_LIST.dat', dtype=test_conf,skiprows=1,ndmin=1)

for t,test in enumerate(TEST_LIST):

  if (t%nranks == rank):

     dba.do_big_ave(test,POSTPROC_DIR)
     dbaHF.do_big_ave_HF(test,POSTPROC_DIR)
     dbap.do_big_ave_phys(test,POSTPROC_DIR)
#    dbad.do_big_ave_dia(test,POSTPROC_DIR)
#    dbaf.do_big_ave_flux(test,POSTPROC_DIR)
#    plot_PTSW.plot_hovmoeller_PTSW(test,POSTPROC_DIR)
#    plot_LDNCCC.plot_hovmoeller_LDNCCC(test,POSTPROC_DIR)
#    plot_LDNC.plot_hovmoeller_LDNC(test,POSTPROC_DIR)
#    plot_CCCCBBBB.plot_hovmoeller_CCCCBBBB(test,POSTPROC_DIR)

     for DATADIR in LIST_POSTPROC:
         plot_C_CF.plot_hovmoeller_C_CF(test,DATADIR,0)
#    plt_BUD.plot_hovmoeller_BUDGET(test,POSTPROC_DIR)
