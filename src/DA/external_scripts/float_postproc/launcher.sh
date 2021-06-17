#! /bin/bash

module purge
module load profile/advanced
module load autoload intelmpi/5.0.1--binary netcdf/4.1.3--intel--cs-xe-2015--binary
module load autoload python/2.7.9 numpy/1.9.2--python--2.7.9 matplotlib/1.4.3--python--2.7.9 scipy/0.15.1--python--2.7.9
# source /pico/work/OGS16_PRACE_P/COPERNICUS/sequence.sh
source /gpfs/work/OGS16_PRACE_P/COPERNICUS/py_env_2.7.9/bin/activate
PYTHONPATH=$PYTHONPATH:/gpfs/work/OGS16_PRACE_P/COPERNICUS/bit.sea

#idate=20150106
#for idate in 20150106 20150113 20150120 20150127
for idate in 20150804 20150811 20150818 20150825
  do

    INPUTDIR=/pico/scratch/userexternal/lmariott/FLOAT_DA_01/wrkdir/POSTPROC/output/AVE_FREQ_1/TMP/
    BASEDIR=/pico/scratch/userexternal/lmariott/FLOAT_DA_01/wrkdir/float_preproc_lm/PROFILATORE_WEEKLY_LOV_OGSTM/
    export MASKFILE=/pico/home/usera07ogs/a07ogs00/OPA/V2C/etc/static-data/MED1672_cut/MASK/meshmask.nc
    DESTDIR=/pico/scratch/userexternal/lmariott/FLOAT_DA_01/wrkdir/MODEL/DA__FREQ_1/

    python float_extractor.py  -t $idate -i $INPUTDIR -b $BASEDIR
#python preproc.py          -t $idate -i $INPUTDIR -b $BASEDIR
#cp $idate.P_l_arg_mis.dat $DESTDIR

    python postproc_float_3dvar.py -t $idate -i $INPUTDIR -b $BASEDIR

    echo " "
done


