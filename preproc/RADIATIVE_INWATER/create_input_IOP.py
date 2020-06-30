#!/bin/env ipython

from basins import V2 as OGS
from checkvar import *
from create_IOP import *
from commons.layer import Layer
from instruments import optbio_float_2019
from instruments import var_conversions
from instruments.matchup_manager import Matchup_Manager
from matchup import statistics
from model_params import *
import numpy as np
import os, sys
from profiler import *
import scipy.io.netcdf as NC
import subprocess
import time
from write_matchup import *

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

''''
Read model data that corresponds to the position and time of BGC-Argo floats
'''

if len(sys.argv) != 3:
    raise ValueError('Wrong inputs!')

M = Matchup_Manager(ALL_PROFILES,TL,BASEDIR)

TI = TimeInterval(sys.argv[1], sys.argv[2],'%Y%m%d-%H:%M:%S')

variable='P_l'

Profilelist_aux=optbio_float_2019.FloatSelector(varname,TI , OGS.med)   # len is no. of profiles

# Adjust time from UTC to local!
for FLOAT in Profilelist_aux:
    FLOAT.time += timedelta(hours=24./360.*FLOAT.lon)

Profilelist_aux2     = [p for p in Profilelist_aux if TI.contains(p.time)]  # additional check
            
Profilelist          = [p for p in Profilelist_aux2 if (int(p.time.strftime('%H'))>10 and int(p.time.strftime('%H'))<14 )]

for p in Profilelist:#[rank::nranks]:
    profile_ID = p.ID()
    print(profile_ID)
    
    bool= findVars(p.available_params)
    if bool == False:
        print('Variables missing for further calculations')
        continue
    ''' 
    phase 1. write OASIM.txt file 
    '''
    
    List_Ed = [M.getMatchups([p], nav_lev, modelvar).subset(Layer(0,1.5)) for modelvar in str_Ed]
    List_Es = [M.getMatchups([p], nav_lev, modelvar).subset(Layer(0,1.5)) for modelvar in str_Es]
    
    Ed = np.asarray([0. if len(List_Ed[i].Model)==0 else List_Ed[i].Model[0] for i in range(len(List_Ed))])
    Es = np.asarray([0. if len(List_Es[i].Model)==0 else List_Es[i].Model[0] for i in range(len(List_Ed))])
    
    if Ed.all() == 0. and Es.all() == 0.:
        print('No model data for this profile')
        continue
    
    if (Ed[4:9].max() + Es[4:9].max()) < 30.:
        print('Low irradiance values of OASIM!')
        continue

    '''
    phase 2. Read BGC-ARGO profiles
    '''
    PresCHL, CHLz,    Qc = p.read('CHL')
    Pres380, Ed_380,  Qc = p.read('IRR_380')
    Pres412, Ed_412,  Qc = p.read('IRR_412')
    Pres490, Ed_490,  Qc = p.read('IRR_490')
    PresPAR, PAR,     Qc = p.read('PAR')
    Lon = p.lon
    Lat = p.lat
    timestr = p.time.strftime("%Y%m%d-%H:%M:%S")
    nLevels = len(PresCHL)
    init_rows = str(timestr) + '\n' + str(Lat) + '\n' + str(nLevels)
    
    if PresCHL[0] == 0.:
        print('First depth equals 0')
        continue
    
    if Ed_380[0] < 30. or Ed_412[0] < 30. or Ed_490[0] < 30:
        print('BGC-Argo low irradiance values - cloud coverage')
        continue
    
    if PresCHL.max() < 15:
        print('Depth range too small')
        continue

    '''
    phase 3. Calculate and save IOPs  
    '''
    PFT1, PFT2, PFT3, PFT4 = PFT_calc(CHLz, 0.40, 0.30, 0.25, 0.05)#0.30, 0.20, 0.40, 0.10)
    #PFT1, PFT2, PFT3, PFT4 = PFT_MED(CHLz)
    
    aNAP  = NAP_abs( CHLz,   0.0129, 0.00862)
    aCDOM = CDOM_abs(CHLz,   0.015, 0.05)
    
    file_cols_PFT = np.vstack((PresCHL, PFT1, PFT2, PFT3, PFT4)).T
    np.savetxt(profile_ID + '_PFT.txt', file_cols_PFT, header = init_rows, delimiter='\t', comments='')
    
    Pres = PresCHL.reshape(PresCHL.shape[0], 1)
    
    file_cols_CDOM = np.hstack((Pres, aCDOM))
    np.savetxt(profile_ID + '_CDOM.txt', file_cols_CDOM, delimiter='\t', comments='' )
    
    file_cols_NAP = np.hstack((Pres, aNAP))
    np.savetxt(profile_ID + '_NAP.txt', file_cols_NAP, delimiter='\t', comments='' )
    
    floatname = profile_ID + '.nc'
    
    np.savetxt(profile_ID + '_OASIM.txt', np.c_[Ed, Es])
    
    '''  
    phase 4 : Run Fortran code
    '''
    command='./compute_IOP.xx ' + profile_ID + '_OASIM.txt ' + profile_ID + '_PFT.txt '  + profile_ID + '_CDOM.txt '  + profile_ID + '_NAP.txt ' + str(floatname) + ' >> log'
    print rank, command 
    subprocess.call(command, shell=True)
    
    '''
    phase 5: Prepare irradiance output .nc files for ARGO-model matchup
    '''    
    ncin=NC4.Dataset(floatname,"r")
    
    Ed380_model  =  np.array( 0.8  * (ncin.variables['Edz'][3,1:] + ncin.variables['Esz'][3,1:])  + 0.2 *  (ncin.variables['Edz'][4,1:] + ncin.variables['Esz'][4,1:]))  * 4 # = 10**(-6) / (10**(-4) * 25) 
    Ed412_model  =  np.array( 0.52 * (ncin.variables['Edz'][4,1:] + ncin.variables['Esz'][4,1:])  + 0.48 * (ncin.variables['Edz'][5,1:] + ncin.variables['Esz'][5,1:]))  * 4 #  W/m2 to muW/cm2
    Ed490_model  =  np.array( 0.4  * (ncin.variables['Edz'][7,1:] + ncin.variables['Esz'][7,1:])  + 0.6 *  (ncin.variables['Edz'][8,1:] + ncin.variables['Esz'][8,1:]))  * 4
    
    ncin.close()
    '''Interpolate Ed380 on CHL (OASIM model) depth quotes'''
    
    Ed380_float = np.interp(PresCHL, Pres380, Ed_380)
    Ed412_float = np.interp(PresCHL, Pres412, Ed_412)
    Ed490_float = np.interp(PresCHL, Pres490, Ed_490)
    
    ncout = save_matchup(floatname, PresCHL, Ed380_float, Ed412_float, Ed490_float, Ed380_model, Ed412_model, Ed490_model, timestr)
    
    '''Move the in-water radiative transfer model output to a separate directory'''
    movefiles = 'mv ' + str(floatname) + ' NCOUT/'
    os.system(movefiles)
    
    
    ''' Move the .txt files you don't need any more '''
    txtfiles1 = 'mv ' + profile_ID + '_OASIM.txt' + ' TXT_FILES/' 
    txtfiles2 = 'mv ' + profile_ID + '_PFT.txt'   + ' TXT_FILES/' 
    txtfiles3 = 'mv ' + profile_ID + '_NAP.txt'   + ' TXT_FILES/' 
    txtfiles4 = 'mv ' + profile_ID + '_CDOM.txt'  + ' TXT_FILES/'
    os.system(txtfiles1)
    os.system(txtfiles2)
    os.system(txtfiles3)
    os.system(txtfiles4)
