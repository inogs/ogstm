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
    os.system("ln -fs  " + CODEPATH +  "OGSTM_BUILD_DBG/ogstm.xx "+ test['Dir'].decode() + "/" )
    
    os.system("cp subgen.py " + test['Dir'].decode() )
    os.system("cp boundaries.nml " + test['Dir'].decode() )
    os.system("cp oasim_config.yaml " + test['Dir'].decode() )
    os.system("cp -r bcs " + test['Dir'].decode() )
    print("edit namelist.init to set :") 
    print("   variable_rdt  = .false.")
    print("   mld_flag = .false.")
#   print("   ingv_files_direct_reading = .false.")
    print("   ingv_lon_shift = 0 ")
    print("   is_free_surface = .false.")

    if test['BGC_TYPE'].decode() in ['DEFAULT','Default','default']:
        os.system("cp -pf namelist.init.default  " +  test['Dir'].decode() + "/namelist.init")
        os.system("cp -pf namelist.passivetrc.default  " +  test['Dir'].decode() + "/namelist.passivetrc")

    elif test['BGC_TYPE'].decode() in ['BFM','bfm']:
        namelists= CODEPATH +  "/ogstm/ready_for_model_namelists/* "
        os.system("cp -pf " + namelists + test['Dir'].decode() + "/")

    elif test['BGC_TYPE'].decode() in ['FABM-BFM','fabm-bfm']:
        os.system("cp -pf namelist.init.fabm  " +  test['Dir'].decode() + "/namelist.init")

        namelist_optics = CODEPATH +  "/ogstm/ready_for_model_namelists/namelist.optics "
        os.system("cp -pf " + namelist_optics +  test['Dir'].decode() + "/")

        namelist_phys   = CODEPATH +  "/ogstm/ready_for_model_namelists/namelist.phys "
        os.system("cp -pf " + namelist_phys +  test['Dir'].decode() + "/")

        fabm_yaml=  CODEPATH + "/fabm/extern/ogs/fabm_multispectral_2xDetritus.yaml "
        os.system("cp -pf " + fabm_yaml + test['Dir'].decode() + "/fabm.yaml")
        
        ret = os.system("cp -pf ogstm.yaml " + test['Dir'].decode() + "/ogstm.yaml")
        if ret != 0:
            print("ogstm.yaml not found, creating template...")
            os.system("python create_ogstm_yaml.py")
            os.system("cp -pf ogstm.yaml " + test['Dir'].decode() + "/ogstm.yaml")

    elif test['BGC_TYPE'].decode() in ['FABM-ROSENMCARTUR','fabm-rosenmcartur']:
        os.system("cp -pf namelist.init.fabm  " +  test['Dir'].decode() + "/namelist.init")

        namelist_optics = CODEPATH +  "/ogstm/ready_for_model_namelists/namelist.optics "
        os.system("cp -pf " + namelist_optics +  test['Dir'].decode() + "/")

        namelist_phys   = CODEPATH +  "/ogstm/ready_for_model_namelists/namelist.phys "
        os.system("cp -pf " + namelist_phys +  test['Dir'].decode() + "/")

        fabm_yaml= "fabm_rosenmcartur.yaml "
        os.system("cp -pf " + fabm_yaml + test['Dir'].decode() + "/fabm.yaml")
        
        ret = os.system("cp -pf ogstm.yaml " + test['Dir'].decode() + "/ogstm.yaml")
        if ret != 0:
            print("ogstm.yaml not found, creating template...")
            os.system("python create_ogstm_yaml.py")
            os.system("cp -pf ogstm.yaml " + test['Dir'].decode() + "/ogstm.yaml")
    else:

        print("BGC_TYPE " + test['BGC_TYPE'].decode() + "wrong/undefined choose within [DEFAULT[Default,default], BFM[bfm], FABM[fabm]]")
        sys.exit()

