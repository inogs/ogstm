#!/bin/env python

from __future__ import print_function
from datetime import timedelta
import numpy as np
import scipy.io.netcdf as NC

from basins.region import Rectangle
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from instruments.matchup_manager import Matchup_Manager
#from instruments import optbio_float_2019 as optbio_float
from static import superfloat as optbio_float


INPUTDIR  = '/gpfs/scratch/userexternal/eterzic0/OASIM_HF_INWATER/AVEDATA/'
BASEDIR   = '/gpfs/scratch/userexternal/eterzic0/OASIM_HF_INWATER/PROFILATORE/'

START_DATE = '20120101-00:00:00'
END___DATE = '20171231-00:00:00'

TI    = TimeInterval(START_DATE, END___DATE, '%Y%m%d-%H:%M:%S')
TL    = TimeList.fromfilenames(TI, INPUTDIR,"ave*.nc",filtervar="Ed_400")


#export MASKFILE=$PWD/meshmask.nc

ALL_PROFILES = optbio_float.FloatSelector(None, TI, Rectangle(-6,36,30,46))

for FLOAT in ALL_PROFILES:
	FLOAT.time += timedelta(hours=24./360.*FLOAT.lon)  # Change from GMT to local time

vardescriptorfile='/galileo/home/userexternal/eterzic0/CODE/ogstm/preproc/RADIATIVE_INWATER/OASIM_postproc/VarDescriptorB.xml'

#This previous part will be imported in matchups setup.

# The following part, the profiler, is executed once and for all.
# It might take some time, depending on length of simulation or size of files.
if __name__ == '__main__':
    # Here instruments time and positions are read as well as model times
    M = Matchup_Manager(ALL_PROFILES,TL,BASEDIR)

    profilerscript = BASEDIR + 'jobProfiler.sh'
    aggregatedir="/gpfs/scratch/userexternal/eterzic0/OASIM_HF__INWATER/AVEDATA/"
    M.writefiles_for_profiling(vardescriptorfile, profilerscript, aggregatedir=aggregatedir) # preparation of data for aveScan

    M.dumpModelProfiles(profilerscript) # sequential launch of aveScan
