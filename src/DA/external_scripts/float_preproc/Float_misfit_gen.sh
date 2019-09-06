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
module load profile/base
#module load intel/pe-xe-2018--binary
module load autoload
#module load netcdf/4.6.1--intel--pe-xe-2018--binary
module load python/2.7.12 #scipy/1.1.0--python--2.7.12
module load intelmpi/2018--binary

source /gpfs/work/OGS18_PRACE_P_0/COPERNICUS/py_env_2.7.12/bin/activate
#source /gpfs/work/OGS18_PRACE_P_0/COPERNICUS/sequence.sh
export PYTHONPATH=$PYTHONPATH:/gpfs/work/OGS18_PRACE_P_0/COPERNICUS/bit.sea

export OPA_HOME=/gpfs/scratch/userexternal/ateruzzi/DA_Float/RUN_FLOAT_chl_n/


      DA_DIR=$OPA_HOME/wrkdir/MODEL/DA__FREQ_1/
     TMP_DIR=$OPA_HOME/wrkdir/POSTPROC/output/DA__FREQ_1/TMP
  CHLSUP_DIR=$OPA_HOME/wrkdir/POSTPROC/output/DA__FREQ_1/CHL_SUP
     BASEDIR=$DA_DIR/PROFILATORE_WEEKLY_LOV_OGSTM/
      OPADIR=$OPA_HOME/wrkdir/float_preproc/
    DEST_DIR=$OPA_HOME/wrkdir/MODEL/DA__FREQ_1/

export ONLINE_REPO=/gpfs/scratch/userexternal/ateruzzi/REPO_LOV/ONLINE/
export MASKFILE=$OPA_HOME/wrkdir/MODEL/meshmask.nc
mkdir -p $BASEDIR
mkdir -p $DA_DIR/links
mv $DA_DIR/RSTbefore*00000* $DA_DIR/links

DADEP=600

for vv in N3n P_l; do

echo --- Variable $vv
rm -rf $BASEDIR
mkdir -p $BASEDIR
cd $DIR

if [ $vv == P_l ]; then

    mkdir -p $TMP_DIR
    mkdir -p $CHLSUP_DIR

    echo --- var_aggregator.py -l RSTbefore.${DATE}*P1l.nc -i $DA_DIR -d VarDescriptor_P_lagg.xml -t $TMP_DIR  -c $CHLSUP_DIR -m $MASKFILE
    python var_aggregator.py -l RSTbefore.${DATE}*00:00*P1l.nc -i $DA_DIR -d VarDescriptor_P_lagg.xml -t $TMP_DIR  -c $CHLSUP_DIR -m $MASKFILE
    MODEL_AVEDIR=$TMP_DIR
    if [ $? -ne 0 ] ; then exit 1 ; fi
fi

if [ $vv == N3n ] ; then
    MODEL_AVEDIR=$DA_DIR
fi

echo model_avedir $MODEL_AVEDIR

echo --- float_extractor.py -t ${DATE}  -i $MODEL_AVEDIR -b $BASEDIR  -d $OPADIR -v $vv
python float_extractor.py -t ${DATE}  -i $MODEL_AVEDIR -b $BASEDIR  -d $OPADIR -v $vv
if [ $? -ne 0 ] ; then exit 1 ; fi

echo --- preproc.py -t ${DATE}  -i $MODEL_AVEDIR -b $BASEDIR -m $MASKFILE -v $vv -d $DADEP
python preproc.py              -t ${DATE}  -i $MODEL_AVEDIR -b $BASEDIR -m $MASKFILE -v $vv -d $DADEP
if [ $? -ne 0 ] ; then exit 1 ; fi

#echo --- mv $DA_DIR/links/* $DA_DIR
echo ---  copy files cp ${DATE}.${vv}_arg_mis.dat $DEST_DIR
cp ${DATE}.${vv}_arg_mis.dat $DEST_DIR
if [ $? -ne 0 ] ; then exit 1 ; fi

done

mv $DA_DIR/links/RSTbefore* $DA_DIR

echo merge arg_mis and copy
echo ---  python merge_arg_mis.py -t ${DATE}
python merge_arg_mis.py -t ${DATE}
cp ${DATE}.arg_mis.dat $DEST_DIR




