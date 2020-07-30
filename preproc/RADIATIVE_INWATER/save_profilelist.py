#!/bin/env python

''' 
WARNING !!! Before you run this code make sure you load the env.sh !!! 

Save profile and floatlists as .pkl files to facilitate the reading for each simulation 
'''

from __future__ import print_function
from datetime import timedelta
import pickle 

from basins import V2 as OGS
#from instruments import optbio_float_2019 as optbio_float
from instruments import var_conversions
from static import superfloat as optbio_float

from ancillary import *
from configuration import *


variable='P_l'  # The variable of reference - in this case Chl, which is present in all profiles

#varname=var_conversions.FLOAT_OPT_VARS_2019[variable]
varname=var_conversions.SUPERFLOAT_VARS[variable]

Profilelist_aux=optbio_float.FloatSelector(varname, TI ,OGS.med)   # len is no. of profiles

print('Length of the initial profile list = ', len(Profilelist_aux))

Profilelist = []

ctime, ctint, cvars = 0, 0, 0

for p in Profilelist_aux:
	#p.time += timedelta(hours=24./360.* p.lon)   # Adjust time from UTC to local!
	if not TI.contains(p.time): ctime +=1 ; continue
	if not ( int(p.time.strftime('%H'))>10 and int(p.time.strftime('%H'))<14 ): ctint +=1; continue
	#if not findVars(p.available_params.split(' ')[1:], allvars=['TEMP', 'SALI', 'CHLA']): cvars +=1; continue
	#if not findVars(p.available_params.split(' ')[1:], allvars=['TEMP', 'SALI', 'CHLA', 'IRR_380', 'IRR_412', 'IRR_490']): cvars +=1; continue
	if not findVars(p.available_params.split(' ')[1:], allvars=['TEMP', 'SALI', 'CHLA', 'IRR_380', 'IRR_412', 'IRR_490', 'BBP700', 'CDOM' ]): cvars +=1; continue

	#if not findVars(p.available_params.split(' ')[1:], allvars=['TEMP', 'SALI', 'CHLA', 'IRR_380', 'IRR_412', 'IRR_490', 'BBP700', 'CDOM']): cvars +=1; continue
	Profilelist.append(p)

print('Counters: ', ctime, ctint, cvars)

print('Length of the filtered profile list = ', len(Profilelist))

Floatlist=optbio_float.get_wmo_list(Profilelist)               # len is no. of floats


output = open('Profilelist.pkl', 'wb')

pickle.dump(Profilelist, output)
pickle.dump(Floatlist,   output)

output.close()

