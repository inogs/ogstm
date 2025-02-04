import argparse
from bitsea.utilities.argparse_types import existing_dir_path, existing_file_path, generic_path
from bitsea.utilities.argparse_types import date_from_str

def argument():
    parser = argparse.ArgumentParser(description='''
    Writes misfit file and check files in OUTDIR
    ''',
                                     formatter_class=argparse.RawTextHelpFormatter
                                     )

    parser.add_argument('--time', '-t',
                        type=date_from_str,
                        required=True,
                        help='Input time in yyyymmdd format')
    parser.add_argument('--inputdir', '-i',
                        type=existing_dir_path,
                        required=True,
                        help='input dir validation')
    parser.add_argument('--maskfile', '-m',
                        type=existing_file_path,
                        required=True)
    parser.add_argument('--basedir', '-b',
                        type=existing_dir_path,
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
                        type=generic_path,
                        default=None,
                        required=True,
                        help=''' output misfit file ''')
    parser.add_argument('--outdir','-o ',
                        type=existing_dir_path,
                        default=None,
                        required=True,
                        help=''' output directory of check files''')
    return parser.parse_args()


args = argument()

from bitsea.instruments.matchup_manager import Matchup_Manager
from bitsea.commons.mask import Mask
import bitsea.basins.OGS as OGS
from bitsea.instruments.var_conversions import FLOATVARS
from datetime import datetime
from bitsea.commons.layer import Layer
from bitsea.commons.Timelist import TimeList
from bitsea.commons import timerequestors
import os,sys
import numpy as np
from bitsea.instruments import check

profilesource=os.getenv("PROFILES_SOURCE")
if profilesource is None:
    print("Error: Environment variables PROFILES_SOURCE - superfloat or ppcon - must be defined.")
    sys.exit(1)
assert profilesource in ["superfloat", "ppcon"]
if profilesource=="superfloat":
    from bitsea.instruments import superfloat as bio_float
if profilesource=="ppcon":
    if args.variable == 'P_l':
        from bitsea.instruments import superfloat as bio_float ## ppcon chl not yet tested for DA
    else:
        from bitsea.instruments import float_ppcon as bio_float


datestr   = args.time.strftime("%Y%m%d")
BASEDIR  = args.basedir
INPUTDIR = args.inputdir
OUTDIR   = args.outdir
varmod   = args.variable
TheMask = Mask(args.maskfile)
deplim  = int(args.deplim)

nav_lev = TheMask.zlevels
layer   = Layer(0, deplim)

errbase = 0.24  # (Mignot et al., 2019)
errbase_ppcon = [0.44, 0.69, 0.61] #0-200 200-400 and 400-600
erro2obase = 5.  # (Approximation based on QuID V7c evaluation)
Check_Obj = check.check(OUTDIR,verboselevel=1,threshold_nitrate=2)

year  = args.time.year
month = args.time.month
day   = args.time.day

req = timerequestors.Daily_req(year, month, day)
TI = req.time_interval
TI.end_time = datetime(year, month, day, 23, 59)

errorfloat = [0.0690, 0.0969, 0.0997, 0.0826, 0.0660, 0.0500, 0.0360, 0.0140, 0.0320, 0.0390, 0.0340, 0.0490]

#errorfloat_ppcon = [ 0.12  , 0.14  , 0.14  , 0.84  , 0.67  , 0.14  , 0.08  , 0.07  , 0.08  , 0.06  , 0.06  , 0.07  ] # future dev
# method: desroziers 2005 [err float + extra uncert. ]
DICTflag = {'P_l': 0, 'N3n': 1, 'O2o': 2}

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
        print("High frequency float: %d profiles on same day for wmo %s " %(nTimes, p._my_float.wmo))
        TimeRef=datetime(p.time.year,p.time.month,p.time.day,13)
        times= [pp.time for pp in float_track_list]
        delta_t = [t-TimeRef for t in times]
        deltadays = np.zeros((nTimes,) , np.float32)
        for i,d in enumerate(delta_t):
            deltadays[i] = d.days + d.seconds/(86400.)
        chosen=np.abs(deltadays).argmin()
        print("chosen: ", times[chosen])
        p = float_track_list[chosen]
        for pp in float_track_list:
            if not pp==p:
                print("Rejected : ", p.ID(), pp.time)
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
                Goodlist.append(p)


    print("number of Goodlist Profiles:",  len(Goodlist))
    if (Goodlist != []):
        if np.any(LISTexclusion):
           mindepexc = np.min(LISTdepthexc)
           deplim = np.min([mindepexc,deplim])
        ii = Pres <= deplim
        nLevels = ii.sum()
        levels = Pres[ii]
        sfm = singlefloatmatchup.subset(Layer(0,deplim))
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
                if (Qc <0).any():
                   if   lev <=200: errnut = errbase_ppcon[0]
                   elif (lev >=400) and (lev <450): errnut = errbase_ppcon[2] # layer 400m-450m no localization
                   elif (lev >=450) and (lev <600): errnut = errbase_ppcon[2]*(1+(1.5-1.)/(600-450)*(lev-450)) #localization
                   else:  errnut = errbase_ppcon[1]
                else:
                   if lev <= 450: errnut = errbase
                   if (lev > 450):  errnut = errbase*(1+(1.5-1.)/(600-450)*(lev-450))
                errnut = errbase
                testo += "\t%10.5f\t%i\n" % (errnut, Profile[ ilev, 5])
            if varmod == 'O2o':
                erro2o = erro2obase
                testo += "\t%10.5f\t%i\n" % (erro2o, Profile[ ilev, 5])
            if varmod == 'P_l':
                testo += "\t%10.5f\t%i\n" % (errorfloat[month-1], Profile[ilev, 5])
            MISFIT_LINES.append(testo)

with open(args.misfit,'w') as f:
    topstr = "%i\n" % len(MISFIT_LINES)
    f.write(topstr)
    f.writelines(MISFIT_LINES)

checkfile_txt = OUTDIR / (datestr + varmod + '_check.txt')
with open(checkfile_txt,'wt') as fid:
    fid.writelines(LINES)

