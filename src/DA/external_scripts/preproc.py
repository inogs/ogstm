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
from commons.Timelist import TimeList
from commons.layer import Layer
import pylab as pl
from commons.utils import addsep
import basins.OGS as OGS
from instruments import lovbio_float as bio_float
from instruments.var_conversions import LOVFLOATVARS
from commons.mask import Mask
from instruments.matchup_manager import Matchup_Manager
import sys
#import makeplots
import qualitycheck


idate0=args.inputdate
BASEDIR=addsep(args.basedir)
INPUTDIR=addsep(args.inputdir)


# DATE INTERVAL 
year=int(idate0[0:4])
month=int(idate0[4:6])
day=int(idate0[6:8])
idate1=timerequestors.Weekly_req(year,month,day)
print idate1

idate2=timerequestors.Daily_req(year,month,day)
print idate2



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

f = open(idate0+"."+VARLIST[0]+"_arg_mis.dat","w")        # OUTPUT x il 3DVAR

iniz="     \n"
f.writelines(iniz)

# LIST of profiles for the selected variable in VARLIST
# in the interval idate1.time_interval in the mediterranean area 
Profilelist_1=bio_float.FloatSelector(LOVFLOATVARS[VARLIST[0]],idate1.time_interval,OGS.med)
TL = TimeList.fromfilenames(idate2.time_interval, INPUTDIR,"ave.*.nc",filtervar="P_l")
TL.inputFrequency = 'weekly'
M = Matchup_Manager(Profilelist_1, TL ,BASEDIR)

#import sys
#sys.exit()

# LIST of FLOATS (WMO) that are in Profilelist_1
WMOlist=bio_float.get_wmo_list(Profilelist_1)
nWMO=len(WMOlist)              #NUMBER OF float IN THE SELECTED PERIOD

totallines=0
for wmo in WMOlist :
#for ind in np.arange(0,1):
#        wmo=WMOlist[ind]
        print wmo

# SAVE in the "Goodlist" all the profiles for a given FLOAT ID
        Goodlist = []
	SubProfilelist_1 = []
        SubProfilelist_1 = bio_float.filter_by_wmo(Profilelist_1,wmo)
        for i in SubProfilelist_1:
	    Pres, Profile, Qc=i.read('CHLA',read_adjusted[0])   #Profile.shape,Profile.size, np.mean(Profile)
            if(Profile.size!=0) : Goodlist.append(i)     

#  AllProfiles matrix: 0 Model, 1 Ref, 2 Misfit,3 Lat,4 Lon,5 ID FLOAT
#  NOTE: the number of profile per float is len(Goodlist)
        print "number of Goodlist Profiles:",  len(Goodlist)
        #AllProfiles = np.zeros((len(Goodlist),dimnewpress,6),np.float64)   
        OneProfile = np.zeros((1,dimnewpress,6),np.float64)   

        # PROFILES FOR FLOAT AND MODEL
        for ip, pp in enumerate(Goodlist):    
              #print ip,pp.time
              singlefloatmatchup = M.getMatchups([pp], nav_lev, VARLIST[0])
              s200 = singlefloatmatchup.subset(layer)               #values NOT interpolated
              #vertical interpolation at the grid defined in NewPress (1meter)
              s200intmodel=np.interp(NewPres,s200.Depth,s200.Model)  #MODEL (201 values)
              s200intobs_noqc=np.interp(NewPres,s200.Depth,s200.Ref)      #ARGO  (201 values)
              s200intobs=qualitycheck.test(s200intobs_noqc,wmo, pp.lat ,pp.lon)
  
              # TEMPORARY MATRIX WITH ALL THE MATCHUP 
	      if (ip==0):
                 ind=0
                 count=1
                 ix0,jy0=TheMask.convert_lon_lat_to_indices(pp.lon,pp.lat)
                 AllProfiles = OneProfile
                 AllProfiles[ind,:,0] = s200intmodel                     #ONLY MODEL
                 AllProfiles[ind,:,1] = s200intobs                       #ONLY ARGO
                 AllProfiles[ind,:,2] = s200intobs-s200intmodel        #ONLY MISFIT ARGO-MODELLO
                 AllProfiles[ind,:,3] = pp.lat
                 AllProfiles[ind,:,4] = pp.lon
                 AllProfiles[ind,:,5] = wmo                              #ID float
              else:
                 ix1,jy1=TheMask.convert_lon_lat_to_indices(pp.lon,pp.lat)
                 dx=abs(ix1-ix0)
                 dy=abs(jy1-jy0)
                 if (dx<=1. and dy<=1.):
                    count+=1
                    AllProfiles[ind,:,1]+= s200intobs                     #ONLY ARGO
                 else:
                    if (count>1):
                      AllProfiles[ind,:,1]=AllProfiles[ind,:,1]/count                         #ONLY ARGO
                      AllProfiles[ind,:,2]=AllProfiles[ind,:,1]-AllProfiles[ind,:,0]

		    AllProfiles = np.concatenate((AllProfiles,OneProfile),axis=0)
                    ix0,jy0=TheMask.convert_lon_lat_to_indices(pp.lon,pp.lat)
                    count=1
                    ind+=1
                    AllProfiles[ind,:,0] = s200intmodel                     #ONLY MODEL
                    AllProfiles[ind,:,1] = s200intobs                       #ONLY ARGO
                    AllProfiles[ind,:,2] = s200intobs-s200intmodel        #ONLY MISFIT ARGO-MODELLO
                    AllProfiles[ind,:,3] = pp.lat
                    AllProfiles[ind,:,4] = pp.lon
                    AllProfiles[ind,:,5] = wmo                              #ID float

        if (count>1):
            AllProfiles[ind,:,1]=AllProfiles[ind,:,1]/count                         #ONLY ARGO
            AllProfiles[ind,:,2]=AllProfiles[ind,:,1]-AllProfiles[ind,:,0]
      
        totallines+=(AllProfiles.shape[0]*200/5)  
        for jp in  np.arange(0,AllProfiles.shape[0]):
               for lev in range(0,200,5):
                   testo ="\t%i\t%i\t%10.5f\t%10.5f\t%10.5f" %(1,1,AllProfiles[jp,lev,4],AllProfiles[jp,lev,3],lev)
                   testo +="\t%10.5f\t%10.5f" %(0.2, AllProfiles[jp,lev,2] )
                   testo +="\t%10.5f\t%i\n"  %(np.std(AllProfiles[jp,:,2],0),AllProfiles[jp,lev,5])
                   f.writelines(testo)

        #makeplots.plot_floatvsmodel(VARLIST[0],idate0,AllProfiles, NewPres,Goodlist,wmo)
f.close

f = open(idate0+"."+VARLIST[0]+"_arg_mis.dat","w")        # OUTPUT x il 3DVAR
f.seek(0)
iniz="%i" %(totallines)
f.writelines(iniz)
#f.close





