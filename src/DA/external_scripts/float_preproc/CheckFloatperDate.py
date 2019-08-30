# verifica che ci siano i profili dei float e salva la data in daTimes
import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    read something
    ''',
    formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(   '--inputdate', '-t',
                            type = str,
                            required = True,
                            default = 'export_data_ScMYValidation_plan.pkl',
                            help = 'Input date')
    parser.add_argument(   '--inputdir','-i',
                                type = str,
                                required = True,
                                help = 'input dir validation tmp')
    parser.add_argument(   '--basedir', '-b',
                                type = str,
                                default = None,
                                required = True,
                                help = '''output directory, where aveScan.py will run.
                                           Usually called PROFILATORE''')
    return parser.parse_args()

args = argument()

import scipy.io.netcdf as NC
import numpy as np
import os
from commons import timerequestors
from commons.Timelist import TimeList,TimeInterval
from commons.layer import Layer
import pylab as pl
from commons.utils import addsep
import basins.OGS as OGS
from instruments import lovbio_float as bio_float
from instruments.var_conversions import LOVFLOATVARS
from commons.mask import Mask
from instruments.matchup_manager import Matchup_Manager
import sys
import datetime
#import makeplots
import qualitycheck


idate0=args.inputdate
INPUTDIR=addsep(args.inputdir)
BASEDIR=addsep(args.basedir)


# DATE INTERVAL 
#year=int(idate0[0:4])
#month=int(idate0[5:7])
#day=int(idate0[8:10])

year=int(idate0[0:4])
month=int(idate0[4:6])
day=int(idate0[6:8])
#idate1 data del float, cambio per DA ogni 3day
idate1=timerequestors.Interval_req(year,month,day, days=3)
#print idate1.string

#idate2 data del modello
idate2=timerequestors.Daily_req(year,month,day)
#print idate2


# Variable name
VARLIST = ['P_l']      #'N3n','O2o']
read_adjusted = [True] #,False,False]

# MASK of the domain
TheMask=Mask("/pico/home/usera07ogs/a07ogs00/OPA/V2C/etc/static-data/MED1672_cut/MASK/meshmask.nc")
nav_lev = TheMask.zlevels

layer=Layer(0,200)     #layer of the Float profile???? 
# new depth from 0 to 200 meters with 1 meter of resolution 
NewPres=np.linspace(0,200,201)
dimnewpress=len(NewPres)              # data interpolated on the vertical Z

# LIST of profiles for the selected variable in VARLIST
# in the interval idate1.time_interval in the mediterranean area 
Profilelist_1=bio_float.FloatSelector(LOVFLOATVARS[VARLIST[0]],idate1.time_interval,OGS.med)
TL = TimeList.fromfilenames(idate2.time_interval, INPUTDIR,"RST.*.nc",filtervar="P_l",prefix="RST.")
TL.inputFrequency = 'weekly'
M = Matchup_Manager(Profilelist_1, TL ,BASEDIR)


# LIST of FLOATS (WMO) that are in Profilelist_1
WMOlist=bio_float.get_wmo_list(Profilelist_1)
nWMO=len(WMOlist)              #NUMBER OF float IN THE SELECTED PERIOD

for wmo in WMOlist :
        #print wmo

# SAVE in the "Goodlist" all the profiles for a given FLOAT ID
        Goodlist = []
	SubProfilelist_1 = []
        SubProfilelist_1 = bio_float.filter_by_wmo(Profilelist_1,wmo)
        for i in SubProfilelist_1:
            Pres, Profile, Qc=i.read('CHLA',read_adjusted[0])
            if((Pres<200).sum() > 5) :
                      print idate1.string+"-00:00:00"
                      import sys
                      sys.exit()


