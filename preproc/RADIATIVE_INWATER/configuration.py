#!/bin/env python

''' Here you define all the directories, time intervals '''

from __future__ import print_function

from basins.region import Rectangle
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from instruments import optbio_float_2019 as optbio_float


INPUTDIR  = '/gpfs/scratch/userexternal/eterzic0/OASIM_HF_INWATER/AVEDATA/'
BASEDIR   = '/gpfs/scratch/userexternal/eterzic0/OASIM_HF_INWATER/PROFILATORE/'

START_DATE = '20120101-00:00:00'
END___DATE = '20180101-00:00:00'

TI    = TimeInterval(START_DATE, END___DATE, '%Y%m%d-%H:%M:%S')
TL    = TimeList.fromfilenames(TI, INPUTDIR,"ave*.nc",filtervar="Ed_400")

ALL_PROFILES = optbio_float.FloatSelector(None, TI, Rectangle(-6,36,30,46))