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
    parser.add_argument(   '--descriptor', '-d',
                                type = str,
                                default = None,
                                required = True,
                                help = '''VarDescriptor_1.xml, or the complete path''')

    parser.add_argument(   '--variable', '-v',
                                type = str,
                                default = None,
                                required = True,
                                help = '''model variable''')
    return parser.parse_args()


args = argument()

from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from instruments.matchup_manager import Matchup_Manager
import basins.OGS as OGS
from instruments.var_conversions import FLOATVARS
from commons.utils import addsep
from datetime import timedelta
from datetime import datetime
from commons import timerequestors
import numpy as np
import os

profilesource=os.getenv("PROFILES_SOURCE")
if profilesource is None:
    print("Error: Environment variables PROFILES_SOURCE - superfloat or ppcon - must be defined.")
    sys.exit(1)
assert profilesource in ["superfloat", "ppcon"]
if profilesource=="superfloat":
    from instruments import superfloat as bio_float
if profilesource=="ppcon":
    from instruments import float_ppcon as bio_float


meantime=datetime.strptime(args.time,'%Y%m%d')

INPUTDIR=addsep(args.inputdir)
BASEDIR=addsep(args.basedir)


varmod = args.variable

TI_3    = timerequestors.Interval_req(meantime.year, meantime.month, meantime.day, days=1).time_interval
TI_3.end_time = datetime(meantime.year, meantime.month, meantime.day,23,59)

TI_day  = timerequestors.Daily_req( meantime.year, meantime.month, meantime.day).time_interval
TI_day.end_time = datetime(meantime.year, meantime.month, meantime.day,23,59)


Profilelist=bio_float.FloatSelector(FLOATVARS[varmod],TI_3, OGS.med)
print('...              len Profilelist  ' + np.str(len(Profilelist)))
TL = TimeList.fromfilenames(TI_day, INPUTDIR,"RSTbefore*13:00*.nc",filtervar=varmod,prefix="RSTbefore.")
print('...              len TL   ' + np.str(TL.nTimes))
TL.inputFrequency = 'weekly'
M = Matchup_Manager(Profilelist,TL,BASEDIR)
print('...              M=Matchup_Manager done')

profilerscript = BASEDIR + 'jobProfiler.sh'
M.writefiles_for_profiling(args.descriptor, profilerscript, aggregatedir=INPUTDIR) # preparation of data for aveScan
print('...              writesfiles done')

M.dumpModelProfiles(profilerscript) # sequential launch of aveScan
print('...              writesfiles done')



