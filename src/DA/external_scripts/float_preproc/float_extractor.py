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
    parser.add_argument(   '--opadir', '-d',
                                type = str,
                                default = None,
                                required = True,
                                help = '''where I'm working now''')


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
OPADIR =addsep(args.opadir)

#TI_week = timerequestors.Weekly_req(meantime.year, meantime.month, meantime.day).time_interval
# cambio per DA every 3 days
TI_3    = timerequestors.Interval_req(meantime.year, meantime.month, meantime.day, days=1).time_interval
TI_day  = timerequestors.Daily_req( meantime.year, meantime.month, meantime.day).time_interval


Profilelist=bio_float.FloatSelector(None,TI_3, OGS.med)
TL = TimeList.fromfilenames(TI_day, INPUTDIR,"RST*00:00*.nc",filtervar="N1p",prefix="RST.")
TL.inputFrequency = 'weekly'
M = Matchup_Manager(Profilelist,TL,BASEDIR)

profilerscript = BASEDIR + 'jobProfiler.sh'
descriptor=OPADIR+"VarDescriptor.xml"
M.writefiles_for_profiling(descriptor, profilerscript, aggregatedir=INPUTDIR) # preparation of data for aveScan
M.dumpModelProfiles(profilerscript) # sequential launch of aveScan


