import argparse

def argument():
    parser = argparse.ArgumentParser(description='''
    Writes misfit file and  check files in OUTDIR
    ''',
                                     formatter_class=argparse.RawTextHelpFormatter
                                     )

    parser.add_argument('--time', '-t',
                        type=str,
                        required=True,
                        help='Input time in yyyymmdd format')
    parser.add_argument('--inputdir', '-i',
                        type=str,
                        required=True,
                        help='input dir validation')
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
    parser.add_argument('--misfit',
                        type=str,
                        default=None,
                        required=True,
                        help=''' output misfit file ''')
    parser.add_argument('--outdir',
                        type=str,
                        default=None,
                        required=True,
                        help=''' output directory of check files''')
    return parser.parse_args()


args = argument()

from instruments.matchup_manager import Matchup_Manager
from commons.mask import Mask
from instruments import superfloat as bio_float
import basins.OGS as OGS
from instruments.var_conversions import FLOATVARS
from commons.utils import addsep
import datetime
from commons.layer import Layer
from commons.Timelist import TimeList
from commons import timerequestors
import os
import numpy as np
import scipy.io.netcdf as NC
from instruments import check

datestr   = args.time
BASEDIR  = addsep(args.basedir)
INPUTDIR = addsep(args.inputdir)
OUTDIR   = addsep(args.outdir)
varmod   = args.variable
TheMask = Mask(args.maskfile)
deplim  = np.int(args.deplim)

nav_lev = TheMask.zlevels
layer   = Layer(0, deplim)

errbase = 0.24  # (Mignot et al., 2019)
Check_Obj = check.check(OUTDIR,verboselevel=1)

year  = int(datestr[0:4])
month = int(datestr[4:6])
day   = int(datestr[6:8])

req = timerequestors.Daily_req(year, month, day)
TI = req.time_interval
TI.end_time = datetime.datetime(year, month, day, 23, 59)

errorfloat = [0.0690, 0.0969, 0.0997, 0.0826, 0.0660, 0.0500, 0.0360, 0.0140, 0.0320, 0.0390, 0.0340, 0.0490]

DICTflag = {'P_l': 0, 'N3n': 1}

Profilelist = bio_float.FloatSelector(FLOATVARS[varmod], TI, OGS.med)
TL = TimeList.fromfilenames(TI, INPUTDIR,"RSTbefore.*13:00*.nc", filtervar=varmod, prefix="RSTbefore.")
TL.inputFrequency = 'daily'
M = Matchup_Manager(Profilelist, TL, BASEDIR)

WMOlist = bio_float.get_wmo_list(Profilelist)

LINES = []
MISFIT_LINES=[]
for wmo in WMOlist:
    Goodlist      = []# SAVE in the "Goodlist" all the profiles for a given FLOAT ID
    LISTexclusion = []
    LISTdepthexc  = []
    float_track_list = bio_float.filter_by_wmo(Profilelist, wmo)

    for i in float_track_list:
        Pres, Profile, Qc = i.read(FLOATVARS[varmod])
        if((Pres < 200).sum() <= 5): continue

        sfm = M.getMatchups2([i],nav_lev, varmod, checkobj = Check_Obj,interpolation_on_Float=True )
        CheckReport = sfm.CheckReports[0]
        if CheckReport.linestr == '':
            Goodlist.append(i)
        else:
            LINES.append(CheckReport.linestr)
            if varmod == "N3n":
                if CheckReport.reason ==2:
                    LISTexclusion.append(True)
                    LISTdepthexc.append(CheckReport.depthexc)

#  AllProfiles matrix: 0 Model, 1 Ref, 2 Misfit,3 Lat,4 Lon,5 ID FLOAT
#  NOTE: the number of profile per float is len(Goodlist)
    print "number of Goodlist Profiles:",  len(Goodlist)
    if (Goodlist != []):
        if np.any(LISTexclusion):
           mindepexc = np.min(LISTdepthexc)
           deplim = np.min([mindepexc,deplim])
        ii = Pres < deplim
        nLevels = ii.sum()
        OneProfile = np.zeros((1, nLevels, 6), np.float64)

        AllProfiles = OneProfile # TEMPORARY MATRIX WITH ALL THE MATCHUP
        levels = Pres[ii]
        print ' levels shape ' , nLevels
        for ip, pp in enumerate(Goodlist):
            singlefloatmatchup = M.getMatchups([pp], nav_lev, varmod)
            s200 = singlefloatmatchup.subset(layer)  # values NOT interpolated
            s200intmodel = np.interp(levels, s200.Depth, s200.Model)
            s200intobs = Profile[Pres < deplim]

            if (ip == 0):
                count = 1
                ix0, jy0 = TheMask.convert_lon_lat_to_indices(pp.lon, pp.lat)
                AllProfiles[0, :, 0] = s200intmodel  # ONLY MODEL
                AllProfiles[0, :, 1] = s200intobs  # ONLY ARGO
                AllProfiles[0, :, 2] = s200intobs - s200intmodel  # ONLY MISFIT ARGO-MODELLO
                AllProfiles[0, :, 3] = pp.lat
                AllProfiles[0, :, 4] = pp.lon
                AllProfiles[0, :, 5] = wmo
            else:
                ix1, jy1 = TheMask.convert_lon_lat_to_indices(pp.lon, pp.lat)
                dx = abs(ix1-ix0)
                dy = abs(jy1-jy0)
                if (dx <= 1. and dy <= 1.):
                    count += 1
                    AllProfiles[ind, :, 1] += s200intobs  # ONLY ARGO
                else:
                    if (count > 1):
                        AllProfiles[ind, :, 1] = AllProfiles[ind,:, 1]/count  # ONLY ARGO
                        AllProfiles[ind, :, 2] = AllProfiles[ind,:, 1]-AllProfiles[ind, :, 0]

                    AllProfiles = np.concatenate((AllProfiles, OneProfile), axis=0)
                    ix0, jy0 = TheMask.convert_lon_lat_to_indices(pp.lon, pp.lat)
                    count = 1
                    ind += 1
                    AllProfiles[ind, :, 0] = s200intmodel  # ONLY MODEL
                    AllProfiles[ind, :, 1] = s200intobs  # ONLY ARGO
                    AllProfiles[ind, :, 2] = s200intobs - s200intmodel  # ONLY MISFIT ARGO-MODELLO
                    AllProfiles[ind, :, 3] = pp.lat
                    AllProfiles[ind, :, 4] = pp.lon
                    AllProfiles[ind, :, 5] = wmo

        if (count > 1):
            AllProfiles[ind, :, 1] = AllProfiles[ind, :, 1]/count  # ONLY ARGO
            AllProfiles[ind, :, 2] = AllProfiles[ind, :, 1] - AllProfiles[ind, :, 0]

        nProfiles, nLevels, _ = AllProfiles.shape
        if (nProfiles > 0):
            print '   totallines ', nLevels * nProfiles
            for jp in np.arange(nProfiles):
                for ilev in range(nLevels):
                    lev = levels[ilev]
                    testo = "\t%i\t%i\t%10.5f\t%10.5f\t%10.5f" % (1, DICTflag[varmod], AllProfiles[jp, ilev, 4], AllProfiles[jp, ilev, 3], lev)
                    testo += "\t%10.5f\t%10.5f" % (0.2, AllProfiles[jp, ilev, 2])
# errore fisso che aumenta  verso il fondo
                    if varmod == 'N3n':
                        if lev <= 450: errnut = errbase
                        if lev > 450:  errnut = errbase*(1+(1.5-1.)/(600-450)*(lev-450))
                        testo += "\t%10.5f\t%i\n" % (errnut, AllProfiles[jp, ilev, 5])
                    if varmod == 'P_l':
                        testo += "\t%10.5f\t%i\n" % (errorfloat[month-1], AllProfiles[jp, ilev, 5])
                    MISFIT_LINES.append(testo)
        del AllProfiles
        del OneProfile
    else:
        print "Goodlist vuota"
    del Goodlist


f=open(args.misfit,'w')
topstr = "%i\n" % len(MISFIT_LINES)
f.write(topstr)
f.writelines(MISFIT_LINES)
f.close()

checkfile_txt = OUTDIR + time + varmod + '_check.txt'
fid=open(checkfile_txt,'wt')
fid.writelines(LINES)
fid.close()
