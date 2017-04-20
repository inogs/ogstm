import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Creates ave files for aveScan profiler in chain validation
    ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--time','-t',
                                type = str,
                                required = True,
                                help = 'mean time')

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
import matplotlib
matplotlib.use('Agg')
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from instruments.matchup_manager import Matchup_Manager
import basins.OGS as OGS
from instruments import lovbio_float as bio_float
from commons.mask import Mask
from commons.utils import addsep
from datetime import timedelta
from datetime import datetime
from commons import timerequestors


meantime=datetime.strptime(args.time,'%Y%m%d')


INPUTDIR=addsep(args.inputdir)
BASEDIR=addsep(args.basedir)

TI_week = timerequestors.Weekly_req(meantime.year, meantime.month, meantime.day).time_interval
TI_day  = timerequestors.Daily_req( meantime.year, meantime.month, meantime.day).time_interval



Profilelist=bio_float.FloatSelector(None,TI_week, OGS.med)
TL = TimeList.fromfilenames(TI_day, INPUTDIR,"ave*.nc",filtervar="P_l")
TL.inputFrequency = 'weekly'
M = Matchup_Manager(Profilelist,TL,BASEDIR)

profilerscript = BASEDIR + 'jobProfiler.sh'
descriptor="/pico/scratch/userexternal/lmariott/FLOAT_DA_01/wrkdir/float_preproc/VarDescriptor.xml"
M.writefiles_for_profiling(descriptor, profilerscript, aggregatedir=INPUTDIR) # preparation of data for aveScan
M.dumpModelProfiles(profilerscript) # sequential launch of aveScan


