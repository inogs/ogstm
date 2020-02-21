import qualitycheck
import sys
from instruments.matchup_manager import Matchup_Manager
from commons.mask import Mask
from instruments import superfloat as bio_float
import basins.OGS as OGS
from instruments.var_conversions import FLOATVARS
from commons.utils import addsep
import pylab as pl
import datetime
from commons.layer import Layer
from commons.Timelist import TimeList
from commons import timerequestors
import os
import numpy as np
import scipy.io.netcdf as NC
import argparse

import check


def argument():
    parser = argparse.ArgumentParser(description='''
    read something
    ''',
                                     formatter_class=argparse.RawTextHelpFormatter
                                     )

    parser.add_argument('--inputdate', '-t',
                        type=str,
                        required=True,
                        default='export_data_ScMYValidation_plan.pkl',
                        help='Input date')
    parser.add_argument('--inputdir', '-i',
                        type=str,
                        required=True,
                        help='input dir validation tmp')
    parser.add_argument('--maskfile', '-m',
                        type=str,
                        required=True)
    parser.add_argument('--basedir', '-b',
                        type=str,
                        default=None,
                        required=True,
                        help='''output directory, where aveScan.py will run.
                                           Usually called PROFILATORE''')
    parser.add_argument('--variable', '-v',
                        type=str,
                        default=None,
                        required=True,
                        help=''' Model variable''')
    parser.add_argument('--deplim', '-d',
                        type=str,
                        default=None,
                        required=True,
                        help=''' Depth of assimilation''')
    return parser.parse_args()


args = argument()

# import makeplots


idate0 = args.inputdate
BASEDIR = addsep(args.basedir)
INPUTDIR = addsep(args.inputdir)
varmod = args.variable

errbase = 0.24  # (Mignot et al., 2019)

# DATE INTERVAL
year = int(idate0[0:4])
month = int(idate0[4:6])
day = int(idate0[6:8])
# idate1 data del float, cambio per DA ogni 3days
idate1 = timerequestors.Interval_req(year, month, day, days=1)
idate1.time_interval.end_time = datetime.datetime(year, month, day, 23, 59)
print idate1.string

# idate2 data del modello
idate2 = timerequestors.Daily_req(year, month, day)
idate2.time_interval.end_time = datetime.datetime(year, month, day, 23, 59)
print idate2

errorfloat = [0.0690, 0.0969, 0.0997, 0.0826, 0.0660,
              0.0500, 0.0360, 0.0140, 0.0320, 0.0390, 0.0340, 0.0490]


# Variable name
VARLIST = [varmod]  # 'N3n','O2o']
# read_adjusted = [True]  # ,False,False]

DICTflag = {'P_l': 0,
            'N3n': 1}


# MASK of the domain
TheMask = Mask(args.maskfile)
nav_lev = TheMask.zlevels

deplim = np.int(args.deplim)
layer = Layer(0, deplim)  # layer of the Float profile????

f = open(idate0+"."+VARLIST[0]+"_arg_mis.dat", "w")        # OUTPUT x il 3DVAR

iniz = "       \n"
f.writelines(iniz)

# LIST of profiles for the selected variable in VARLIST
# in the interval idate1.time_interval in the mediterranean area
Profilelist_1 = bio_float.FloatSelector(
    FLOATVARS[VARLIST[0]], idate1.time_interval, OGS.med)
TL = TimeList.fromfilenames(idate2.time_interval, INPUTDIR,
                            "RSTbefore.*13:00*.nc", filtervar=varmod, prefix="RSTbefore.")
TL.inputFrequency = 'weekly'
M = Matchup_Manager(Profilelist_1, TL, BASEDIR)


# LIST of FLOATS (WMO) that are in Profilelist_1
WMOlist = bio_float.get_wmo_list(Profilelist_1)
nWMO = len(WMOlist)  # NUMBER OF float IN THE SELECTED PERIOD

print "nWMO", nWMO
totallines = 0
LINES = []
for wmo in WMOlist:
    print wmo

# SAVE in the "Goodlist" all the profiles for a given FLOAT ID
    Goodlist = []
    SubProfilelist_1 = []
    SubProfilelist_1 = bio_float.filter_by_wmo(Profilelist_1, wmo)
    LISTexclusion = []
    LISTdepthexc = []
    for i in SubProfilelist_1:
        # Profile.shape,Profile.size, np.mean(Profile)
        Pres, Profile, Qc = i.read(FLOATVARS[VARLIST[0]])
      #   dimnewpress = Pres[Pres < deplim].shape[0]
        # if(Profile.size!=0) : Goodlist.append(i)
        if((Pres < 200).sum() <= 5): 
           continue
        singlefloatmatchup = M.getMatchups([i], nav_lev, VARLIST[0])
        levelsp = Pres[Pres < deplim]
        s200p = singlefloatmatchup.subset(layer)  # values NOT interpolated
        s200intmodelp = np.interp(levelsp, s200p.Depth, s200p.Model)  # MODEL (201 values)
        s200intobsp = Profile[Pres < deplim]
        if varmod == 'N3n':
           line, depthexc = check.nitrate_check(s200intmodelp, s200intobsp, levelsp, i, './')
           if line != '':
              LINES.append(line)
              reason = np.int(line.rsplit('\t')[-1].rsplit('\n')[0])
              if (reason == 1) or (reason == 3) or (reason == -1):
                 continue
              elif (reason == 2):
                 LISTexclusion.append(True)
                 LISTdepthexc.append(depthexc)
                 Goodlist.append(i)
               #   plotmat[indexexc:, ip] = np.nan
               #   plotmat_model[indexexc:, ip] = np.nan
           else:
              Goodlist.append(i)

        if varmod == 'P_l':
           line = check.chlorophyll_check(s200intmodelp, s200intobsp, levelsp, i, './')
           if line != '':
              LINES.append(line)
              continue
           else:
              Goodlist.append(i)


#  AllProfiles matrix: 0 Model, 1 Ref, 2 Misfit,3 Lat,4 Lon,5 ID FLOAT
#  NOTE: the number of profile per float is len(Goodlist)
    print "number of Goodlist Profiles:",  len(Goodlist)
    if (Goodlist != []):
        if np.any(LISTexclusion):
           mindepexc = np.min(LISTdepthexc)
           deplim = np.min([mindepexc,deplim])
        dimnewpress = Pres[Pres < deplim].shape[0]
        # AllProfiles = np.zeros((len(Goodlist),dimnewpress,6),np.float64)
        OneProfile = np.zeros((1, dimnewpress, 6), np.float64)

        # PROFILES FOR FLOAT AND MODEL
        AllProfiles = OneProfile
        levels = Pres[Pres < deplim]
        print ' levels shape ' + np.str(levels.shape[0])
        for ip, pp in enumerate(Goodlist):
            print ip, pp.time
            singlefloatmatchup = M.getMatchups([pp], nav_lev, VARLIST[0])
            s200 = singlefloatmatchup.subset(layer)  # values NOT interpolated
            s200intmodel = np.interp(
                levels, s200.Depth, s200.Model)  # MODEL (201 values)
            s200intobs = Profile[Pres < deplim]
            # s200intobs = qualitycheck.test(
            #     s200intobs_noqc, wmo, pp.lat, pp.lon)
            

            # TEMPORARY MATRIX WITH ALL THE MATCHUP
            if (ip == 0):
                ind = 0
                count = 1
                ix0, jy0 = TheMask.convert_lon_lat_to_indices(pp.lon, pp.lat)
                AllProfiles[ind, :, 0] = s200intmodel  # ONLY MODEL
                AllProfiles[ind, :, 1] = s200intobs  # ONLY ARGO
                AllProfiles[ind, :, 2] = s200intobs - \
                    s200intmodel  # ONLY MISFIT ARGO-MODELLO
                AllProfiles[ind, :, 3] = pp.lat
                AllProfiles[ind, :, 4] = pp.lon
                AllProfiles[ind, :, 5] = wmo  # ID float
            else:
                ix1, jy1 = TheMask.convert_lon_lat_to_indices(pp.lon, pp.lat)
                dx = abs(ix1-ix0)
                dy = abs(jy1-jy0)
                if (dx <= 1. and dy <= 1.):
                    count += 1
                    AllProfiles[ind, :, 1] += s200intobs  # ONLY ARGO
                else:
                    if (count > 1):
                        AllProfiles[ind, :, 1] = AllProfiles[ind,
                                                             :, 1]/count  # ONLY ARGO
                        AllProfiles[ind, :, 2] = AllProfiles[ind,
                                                             :, 1]-AllProfiles[ind, :, 0]

                    AllProfiles = np.concatenate(
                        (AllProfiles, OneProfile), axis=0)
                    ix0, jy0 = TheMask.convert_lon_lat_to_indices(
                        pp.lon, pp.lat)
                    count = 1
                    ind += 1
                    AllProfiles[ind, :, 0] = s200intmodel  # ONLY MODEL
                    AllProfiles[ind, :, 1] = s200intobs  # ONLY ARGO
                    AllProfiles[ind, :, 2] = s200intobs - \
                        s200intmodel  # ONLY MISFIT ARGO-MODELLO
                    AllProfiles[ind, :, 3] = pp.lat
                    AllProfiles[ind, :, 4] = pp.lon
                    AllProfiles[ind, :, 5] = wmo  # ID float

        if (count > 1):
            AllProfiles[ind, :, 1] = AllProfiles[ind, :, 1]/count  # ONLY ARGO
            AllProfiles[ind, :, 2] = AllProfiles[ind, :, 1] - \
                AllProfiles[ind, :, 0]


        if (AllProfiles.shape[0] > 0):
            totlineswmo = 0
            for iip in range(AllProfiles.shape[0]):
                totlineswmo += AllProfiles[iip, :, :].shape[0]
            # totallines+=(AllProfiles.shape[0]*deplim/5)
            totallines += totlineswmo
            print '   totallines ' + np.str(totallines)
            for jp in np.arange(0, AllProfiles.shape[0]):
                # for lev in range(0,deplim,5):
                for ilev in range(AllProfiles[jp, :, :].shape[0]):
                    lev = levels[ilev]
                    testo = "\t%i\t%i\t%10.5f\t%10.5f\t%10.5f" % (
                        1, DICTflag[VARLIST[0]], AllProfiles[jp, ilev, 4], AllProfiles[jp, ilev, 3], lev)
                    testo += "\t%10.5f\t%10.5f" % (0.2,
                                                   AllProfiles[jp, ilev, 2])
# errore fisso che aumenta  verso il fondo
                    if varmod == 'N3n':
                        if lev <= 450:
                            errnut = errbase
                        if lev > 450:
                            errnut = errbase*(1+(1.5-1.)/(600-450)*(lev-450))
                        testo += "\t%10.5f\t%i\n" % (errnut,
                                                     AllProfiles[jp, ilev, 5])
                    if varmod == 'P_l':
                        testo += "\t%10.5f\t%i\n" % (
                            errorfloat[month-1], AllProfiles[jp, ilev, 5])
                    f.writelines(testo)
        del AllProfiles
        del OneProfile
    else:
        print "Goodlist vuota"
    del Goodlist


f.seek(0)
iniz = "%i" % (totallines)
f.writelines(iniz)
f.close()

fid=open('./OUTTXT/' + idate0 + varmod + '_check.txt','wt')
fid.writelines(LINES)
fid.close()
