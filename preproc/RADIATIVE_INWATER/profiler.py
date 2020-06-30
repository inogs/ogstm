from instruments import optbio_float_2019
from instruments.matchup_manager import Matchup_Manager
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from basins.region import Rectangle

INPUTDIR = '/gpfs/scratch/userexternal/eterzic0/OASIM_HF_INWATER/INDATA/'
BASEDIR='/gpfs/scratch/userexternal/eterzic0/OASIM_HF_INWATER/PROFILATORE/'

DATESTART = '20120101-00:00:00'
DATE__END = '20171231-00:00:00'

#export MASKFILE=$PWD/meshmask.nc

T_INT = TimeInterval(DATESTART,DATE__END, '%Y%m%d-%H:%M:%S')
TL    = TimeList.fromfilenames(T_INT, INPUTDIR,"ave*.nc",filtervar="Ed_380")

ALL_PROFILES = optbio_float_2019.FloatSelector(None,T_INT, Rectangle(-6,36,30,46))

from datetime import timedelta

for FLOAT in ALL_PROFILES:
	FLOAT.time += timedelta(hours=24./360.*FLOAT.lon)  # Change from GMT to local time

vardescriptorfile='/galileo/home/userexternal/eterzic0/CODE/ogstm/preproc/RADIATIVE_INWATER/VarDescriptorB.xml'

#This previous part will be imported in matchups setup.

# The following part, the profiler, is executed once and for all.
# It might take some time, depending on length of simulation or size of files.
if __name__ == '__main__':
    # Here instruments time and positions are read as well as model times
    M = Matchup_Manager(ALL_PROFILES,TL,BASEDIR)

    profilerscript = BASEDIR + 'jobProfiler.sh'
    aggregatedir="/gpfs/scratch/userexternal/eterzic0/OASIM_HF__INWATER/INDATA/"
    M.writefiles_for_profiling(vardescriptorfile, profilerscript, aggregatedir=aggregatedir) # preparation of data for aveScan

    M.dumpModelProfiles(profilerscript) # sequential launch of aveScan
