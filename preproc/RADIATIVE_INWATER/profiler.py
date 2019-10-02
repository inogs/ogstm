from instruments import optbio_float_2019
from instruments.matchup_manager import Matchup_Manager
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from basins.region import Rectangle


INPUTDIR = '/galileo/home/userexternal/eterzic0/OASIM_POSTPROC/ARGO_MATCHUP/INDATA/'
# output directory, where aveScan.py will be run.

BASEDIR='/galileo/home/userexternal/eterzic0/OASIM_POSTPROC/ARGO_MATCHUP/PROFILATORE/'

DATESTART = '20120101'
DATE__END = '20171231'

#export MASKFILE=$PWD/meshmask.nc

T_INT = TimeInterval(DATESTART,DATE__END, '%Y%m%d')
TL = TimeList.fromfilenames(T_INT, INPUTDIR,"ave*.nc",filtervar="Ed_380")

ALL_PROFILES = optbio_float_2019.FloatSelector(None,T_INT, Rectangle(-6,36,30,46))

vardescriptorfile='/galileo/home/userexternal/eterzic0/OASIM_POSTPROC/ARGO_MATCHUP/CODES/VarDescriptorB.xml'

#This previous part will be imported in matchups setup.

# The following part, the profiler, is executed once and for all.
# It might take some time, depending on length of simulation or size of files.
if __name__ == '__main__':
    # Here instruments time and positions are read as well as model times
    M = Matchup_Manager(ALL_PROFILES,TL,BASEDIR)

    profilerscript = BASEDIR + 'jobProfiler.sh'
    aggregatedir="/pico/scratch/userexternal/gbolzon0/eas_v12/eas_v19_3/wrkdir/POSTPROC/output/AVE_FREQ_1/TMP/"
    M.writefiles_for_profiling(vardescriptorfile, profilerscript, aggregatedir=aggregatedir) # preparation of data for aveScan

    M.dumpModelProfiles(profilerscript) # sequential launch of aveScan
