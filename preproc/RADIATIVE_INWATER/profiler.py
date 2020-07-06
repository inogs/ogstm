#!/bin/env python

from __future__ import print_function
from datetime import timedelta

from instruments.matchup_manager import Matchup_Manager

from configuration import *

#export MASKFILE=$PWD/meshmask.nc

ALL_PROFILES = optbio_float.FloatSelector(None, TI, Rectangle(-6,36,30,46))
for FLOAT in ALL_PROFILES:
	FLOAT.time += timedelta(hours=24./360.*FLOAT.lon)  # Change from GMT to local time

vardescriptorfile='/galileo/home/userexternal/eterzic0/CODE/ogstm/preproc/RADIATIVE_INWATER/VarDescriptorB.xml'

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
