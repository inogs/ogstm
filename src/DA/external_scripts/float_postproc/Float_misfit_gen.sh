#! /bin/bash

usage() {
echo "Generates misfit files for argo floats "
echo "SYNOPSYS"
echo "Float_misfit_gen.sh [ -t date]  [ -d dir]"
echo "EXAMPLE"
echo 'Float_misfit_gen.sh -t 20150106 -d .'
echo ""
}


if [ $# -lt 4 ] ; then
  usage
  exit 1
fi

for I in 1 2; do
   case $1 in
      "-t" ) DATE=$2;;
      "-d" ) DIR=$2;;
        *  ) echo "Unrecognized option $1." ; usage;  exit 1;;
   esac
   shift 2
done


module purge
module load profile/advanced
module load autoload intelmpi/5.0.1--binary netcdf/4.1.3--intel--cs-xe-2015--binary
module load autoload python/2.7.9 numpy/1.9.2--python--2.7.9 matplotlib/1.4.3--python--2.7.9 scipy/0.15.1--python--2.7.9
source /gpfs/work/OGS16_PRACE_P/COPERNICUS/py_env_2.7.9/bin/activate
#PYTHONPATH=$PYTHONPATH:/gpfs/work/OGS16_PRACE_P/COPERNICUS/bit.sea
PYTHONPATH=$PYTHONPATH:/pico/scratch/userexternal/lmariott/bit.sea
#PYTHONPATH=$PYTHONPATH:/pico/scratch/userexternal/gbolzon0/bit.sea

export OPA_HOME=FDA_all2015_newstd_3days

##MODEL_AVEDIR=$CINECA_SCRATCH/$OPA_HOME/wrkdir/MODEL/AVE_FREQ_1/
MODEL_AVEDIR=$CINECA_SCRATCH/$OPA_HOME/wrkdir/MODEL/DA__FREQ_1/
     TMP_DIR=$CINECA_SCRATCH/$OPA_HOME/wrkdir/POSTPROC/output/AVE_FREQ_1/TMP
  CHLSUP_DIR=$CINECA_SCRATCH/$OPA_HOME/wrkdir/POSTPROC/output/AVE_FREQ_1/CHL_SUP
     BASEDIR=$CINECA_SCRATCH/$OPA_HOME/wrkdir/float_postproc/PROFILATORE_WEEKLY_LOV_OGSTM/
      OPADIR=$CINECA_SCRATCH/$OPA_HOME/wrkdir/float_postproc/
    DEST_DIR=$CINECA_SCRATCH/$OPA_HOME/wrkdir/MODEL/DA__FREQ_1/
  INDEX_FILE=$CINECA_SCRATCH/$OPA_HOME/wrkdir/float_postproc/$1

     
export MASKFILE=/pico/home/usera07ogs/a07ogs00/OPA/V2C/etc/static-data/MED1672_cut/MASK/meshmask.nc

cd $DIR
#python var_aggregator.py -l RST.${DATE}*P1l.nc -i $MODEL_AVEDIR -d VarDescriptor.xml -t $TMP_DIR  -c $CHLSUP_DIR -m $MASKFILE
if [ $? -ne 0 ] ; then exit 1 ; fi

python float_extractor.py -t ${DATE}  -i $TMP_DIR -b $BASEDIR  -d $OPADIR
if [ $? -ne 0 ] ; then exit 1 ; fi

#python preproc.py              -t ${DATE}  -i $TMP_DIR -b $BASEDIR
python postproc_float_3dvar.py  -t ${DATE}  -i $TMP_DIR -b $BASEDIR -d $DEST_DIR >> $INDEX_FILE
if [ $? -ne 0 ] ; then exit 1 ; fi

#cp ${DATE}.P_l_arg_mis.dat $DEST_DIR
if [ $? -ne 0 ] ; then exit 1 ; fi

