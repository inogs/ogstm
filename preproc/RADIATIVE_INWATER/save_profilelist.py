#!/bin/env python

''' Save profile and floatlists as .pkl files to facilitate the reading for each simulation '''

from __future__ import print_function
from datetime import timedelta
import pickle 

from basins import V2 as OGS
from instruments import optbio_float_2019 as optbio_float
from instruments import var_conversions

from ancillary import *
from configuration import *


variable='P_l'  # The variable of reference - in this case Chl, which is present in all profiles
varname=var_conversions.FLOAT_OPT_VARS_2019[variable]

Profilelist_aux=optbio_float.FloatSelector(varname,TI , OGS.med)   # len is no. of profiles

Profilelist = []

for p in Profilelist_aux:
	p.time += timedelta(hours=24./360.* p.lon)   # Adjust time from UTC to local!
	if not TI.contains(p.time): continue
	if not ( int(p.time.strftime('%H'))>10 and int(p.time.strftime('%H'))<14 ): continue
	if not findVars(p.available_params): continue
	Profilelist.append(p)

Floatlist=optbio_float.get_wmo_list(Profilelist)               # len is no. of floats

output = open('Profilelist.pkl', 'wb')

pickle.dump(Profilelist, output)
pickle.dump(Floatlist,   output)

output.close()

