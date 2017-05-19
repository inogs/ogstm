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
    CODEPATH = test['Code'] + "/ogstm/"
    CODEPATH = CODEPATH.replace("~",os.getenv("HOME"))       
    os.system("ln -fs  " + CODEPATH +  "bin/ogstm.xx "+ test['Dir'] + "/" )

    namelists= CODEPATH +  "ready_for_model_namelists/* "
    os.system("cp -pf " + namelists + test['Dir'] + "/")

