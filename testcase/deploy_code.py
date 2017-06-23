import os,sys
import glob
import numpy as np

from mydtype import *

import scipy.io.netcdf as NC

import pickle

def deploy_code(test):
    os.system("mkdir -p " + test['Dir'].decode())
#######################
#OGSTM
#######################
#   ogstm.xx --> executable
    CODEPATH = test['Code'].decode() 
    CODEPATH = CODEPATH.replace("~",os.getenv("HOME"))
    os.system("ln -fs  " + CODEPATH +  "/OGSTM_BUILD/ogstm.xx "+ test['Dir'].decode() + "/" )

    namelists= CODEPATH +  "/ogstm/ready_for_model_namelists/* "
    os.system("cp -pf " + namelists + test['Dir'].decode() + "/")
    
    os.system("cp subgen.py " + test['Dir'].decode() )
    print("edit namelist.init to set ingv_files_direct_reading = .false. and ingv_lon_shift = 0")
