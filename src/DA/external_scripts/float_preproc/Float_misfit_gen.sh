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
module load autoload
module load intel/pe-xe-2017--binary
module load netcdf/4.4.1--intel--pe-xe-2017--binary
module load python/2.7.12 scipy/0.18.1--python--2.7.12

source /marconi_work/OGS_dev_0/COPERNICUS/py_env_2.7.12/bin/activate
PYTHONPATH=$PYTHONPATH:/marconi_work/OGS_dev_0/COPERNICUS/bit.sea

export OPA_HOME=/marconi_work/OGS_dev_0/MULTIVARIATE/CODE/ogstm/testcase/TEST01


MODEL_AVEDIR=$OPA_HOME/wrkdir/MODEL/DA__FREQ_1/
     TMP_DIR=$OPA_HOME/wrkdir/POSTPROC/output/AVE_FREQ_1/TMP
  CHLSUP_DIR=$OPA_HOME/wrkdir/POSTPROC/output/AVE_FREQ_1/CHL_SUP
     BASEDIR=$OPA_HOME/wrkdir/float_preproc/PROFILATORE_WEEKLY_LOV_OGSTM/
      OPADIR=$OPA_HOME/wrkdir/float_preproc/
    DEST_DIR=$OPA_HOME/wrkdir/MODEL/DA__FREQ_1/

     
export MASKFILE=$OPA_HOME/wrkdir/MODEL/meshmask.nc
mkdir -p $BASEDIR
mkdir -p $MODEL_AVEDIR/links
mv $MODEL_AVEDIR/RST*00000* $MODEL_AVEDIR/links


cd $DIR
#python var_aggregator.py -l RST.${DATE}*P1l.nc -i $MODEL_AVEDIR -d VarDescriptor.xml -t $TMP_DIR  -c $CHLSUP_DIR -m $MASKFILE
#if [ $? -ne 0 ] ; then exit 1 ; fi


python float_extractor.py -t ${DATE}  -i $MODEL_AVEDIR -b $BASEDIR  -d $OPADIR
if [ $? -ne 0 ] ; then exit 1 ; fi

python preproc.py              -t ${DATE}  -i $MODEL_AVEDIR -b $BASEDIR -m $MASKFILE
if [ $? -ne 0 ] ; then exit 1 ; fi

mv $MODEL_AVEDIR/links/* $MODEL_AVEDIR
cp ${DATE}.N3n_arg_mis.dat $DEST_DIR
if [ $? -ne 0 ] ; then exit 1 ; fi

