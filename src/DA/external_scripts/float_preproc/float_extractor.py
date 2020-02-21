import argparse

print '... ... In float_extractor argparse done'

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
    parser.add_argument(   '--opadir', '-d',
                                type = str,
                                default = None,
                                required = True,
                                help = '''where I'm working now''')

    parser.add_argument(   '--variable', '-v',
                                type = str,
                                default = None,
                                required = True,
                                help = '''model variable''')


    return parser.parse_args()

    return parser.parse_args()

args = argument()
#import matplotlib
#matplotlib.use('Agg')
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from instruments.matchup_manager import Matchup_Manager
import basins.OGS as OGS
from instruments import superfloat as bio_float
from instruments.var_conversions import FLOATVARS
#from commons.mask import Mask
from commons.utils import addsep
from datetime import timedelta
from datetime import datetime
from commons import timerequestors
import numpy as np

meantime=datetime.strptime(args.time,'%Y%m%d')
print '... In float_extractor meantime  ' + np.str(meantime)

INPUTDIR=addsep(args.inputdir)
BASEDIR=addsep(args.basedir)
OPADIR =addsep(args.opadir)

varmod = args.variable

#TI_week = timerequestors.Weekly_req(meantime.year, meantime.month, meantime.day).time_interval
# cambio per DA every 3 days
TI_3    = timerequestors.Interval_req(meantime.year, meantime.month, meantime.day, days=1).time_interval
TI_3.end_time = datetime(meantime.year, meantime.month, meantime.day,23,59)

TI_day  = timerequestors.Daily_req( meantime.year, meantime.month, meantime.day).time_interval
TI_day.end_time = datetime(meantime.year, meantime.month, meantime.day,23,59)


Profilelist=bio_float.FloatSelector(FLOATVARS[varmod],TI_3, OGS.med)
print '...              len Profilelist  ' + np.str(len(Profilelist))
TL = TimeList.fromfilenames(TI_day, INPUTDIR,"RSTbefore*13:00*.nc",filtervar=varmod,prefix="RSTbefore.")
print '...              len TL   ' + np.str(TL.nTimes)
TL.inputFrequency = 'weekly'
M = Matchup_Manager(Profilelist,TL,BASEDIR)
print '...              M=Matchup_Manager done'

profilerscript = BASEDIR + 'jobProfiler.sh'
descriptor=OPADIR+"VarDescriptor_" + varmod + ".xml"
M.writefiles_for_profiling(descriptor, profilerscript, aggregatedir=INPUTDIR) # preparation of data for aveScan
print '...              writesfiles done'

M.dumpModelProfiles(profilerscript) # sequential launch of aveScan
print '...              writesfiles done'



