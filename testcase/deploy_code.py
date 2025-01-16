import os,sys
import glob
import numpy as np

from mydtype import *

import scipy.io.netcdf as NC

import pickle

def deploy_code(test):
    os.system("mkdir -p " + test['Dir'])
#######################
#OGSTM
#######################
#   ogstm.xx --> executable
    CODEPATH = test['Code'] 
    CODEPATH = CODEPATH.replace("~",os.getenv("HOME"))
    os.system("ln -fs  " + CODEPATH +  "/OGSTM_BUILD/ogstm.xx "+ test['Dir'] + "/" )

    namelists= CODEPATH +  "/ogstm/ready_for_model_namelists/* "
    os.system("cp -pf " + namelists + test['Dir'] + "/")
    
    os.system("cp subgen.py " + test['Dir'] )
    print("edit namelist.init to set :") 
    print("   ingv_files_direct_reading = .false.")
    print("   ingv_lon_shift = 0 ")
    print("   is_free_surface = .false.")
