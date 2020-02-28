import argparse

def argument():
    parser = argparse.ArgumentParser(description='''
    Writes misfit file and check files in OUTDIR
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
    parser.add_argument('--outdir','-o ',
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

def choose_profile(float_track_list):
    '''
    Chooses one profile from a daily list
    Argument:
    * float_track_list * a list of profile objects for a single day and single wmo
                         It has more than one element when a float has more than one daily cycles.
                         It happens for some high frequency floats that have 2,3,4 daily cycles.
    Returns:
    * p *  a Profile Object, the nearest to assimilation time (13:00)

    '''
    p = float_track_list[0]
    nTimes = len(float_track_list)
    if nTimes > 1 :
        print "High frequency float: %d profiles on same day for wmo %s " %(nTimes, p._my_float.wmo)
        TimeRef=datetime(p.time.year,p.time.month,p.time.day,13)
        times= [pp.time for pp in float_track_list]
        delta_t = [t-TimeRef for t in times]
        deltadays = np.zeros((nTimes,) , np.float32)
        for i,d in enumerate(delta_t):
            deltadays[i] = d.days + d.seconds/(86400.)
        chosen=np.abs(deltadays).argmin()
        print "chosen: ", times[chosen]
        p = float_track_list[chosen]
        for pp in float_track_list:
            if not pp==p:
                print "Rejected : ", p.ID(), pp.time
    return p



LINES = []
MISFIT_LINES=[]
for wmo in WMOlist:
    Goodlist      = [] # SAVE in the "Goodlist" all the profiles for a given FLOAT ID
    LISTexclusion = []
    LISTdepthexc  = []
    float_track_list = bio_float.filter_by_wmo(Profilelist, wmo)
    p = choose_profile(float_track_list)

    Pres, Profile, Qc = p.read(FLOATVARS[varmod])
    if((Pres < 200).sum() <= 5): continue

    singlefloatmatchup = M.getMatchups2([p],nav_lev, varmod, checkobj = Check_Obj,interpolation_on_Float=True )
    CheckReport = singlefloatmatchup.CheckReports[0]
    if CheckReport.line == '':
        Goodlist.append(p)
    else:
        LINES.append(CheckReport.line)
        if varmod == "N3n":
            if CheckReport.reason ==2:
                LISTexclusion.append(True)
                LISTdepthexc.append(CheckReport.depthexc)


    print "number of Goodlist Profiles:",  len(Goodlist)
    if (Goodlist != []):
        if np.any(LISTexclusion):
           mindepexc = np.min(LISTdepthexc)
           deplim = np.min([mindepexc,deplim])
        ii = Pres < deplim
        nLevels = ii.sum()
        levels = Pres[ii]
        sfm = singlefloatmatchup.subset(layer)
        Profile = np.zeros(( nLevels, 6), np.float64)
        Profile[:, 0] = sfm.Model # s200intmodel  # ONLY MODEL
        Profile[:, 1] = sfm.Ref   # s200intobs  # ONLY ARGO
        Profile[:, 2] = sfm.Ref - sfm.Model #s200intobs - s200intmodel  # ONLY MISFIT ARGO-MODELLO
        Profile[:, 3] = p.lat
        Profile[:, 4] = p.lon
        Profile[:, 5] = wmo
        for ilev in range(nLevels):
            lev = levels[ilev]
            testo = "\t%i\t%i\t%10.5f\t%10.5f\t%10.5f" % (1, DICTflag[varmod], Profile[ilev, 4], Profile[ ilev, 3], lev)
            testo += "\t%10.5f\t%10.5f" % (0.2, Profile[ilev, 2])
# errore fisso che aumenta  verso il fondo
            if varmod == 'N3n':
                if lev <= 450: errnut = errbase
                if lev > 450:  errnut = errbase*(1+(1.5-1.)/(600-450)*(lev-450))
                testo += "\t%10.5f\t%i\n" % (errnut, Profile[ ilev, 5])
            if varmod == 'P_l':
                testo += "\t%10.5f\t%i\n" % (errorfloat[month-1], Profile[ilev, 5])
            MISFIT_LINES.append(testo)

f=open(args.misfit,'w')
topstr = "%i\n" % len(MISFIT_LINES)
f.write(topstr)
f.writelines(MISFIT_LINES)
f.close()

checkfile_txt = OUTDIR + datestr + varmod + '_check.txt'
fid=open(checkfile_txt,'wt')
fid.writelines(LINES)
fid.close()
