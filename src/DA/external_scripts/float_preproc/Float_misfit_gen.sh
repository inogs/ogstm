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
module load netcdf/4.6.1--intel--pe-xe-2018--binary
module load numpy/1.15.2--python--2.7.12 
#module load python/2.7.12 #scipy/1.1.0--python--2.7.12
module load intelmpi/2018--binary

source /gpfs/work/OGS20_PRACE_P/COPERNICUS/py_env_2.7.12/bin/activate
#source /gpfs/work/OGS18_PRACE_P_0/COPERNICUS/sequence.sh
export PYTHONPATH=$PYTHONPATH:/gpfs/work/OGS20_PRACE_P/COPERNICUS/bit.sea

export OPA_HOME=/gpfs/scratch/userexternal/ateruzzi/MULTIVARIATE_24/TEST_04/


      DA_DIR=$OPA_HOME/wrkdir/MODEL/DA__FREQ_1/
     TMP_DIR=$OPA_HOME/wrkdir/POSTPROC/output/DA__FREQ_1/TMP
     BASEDIR=$DA_DIR/PROFILATORE_WEEKLY_LOV_OGSTM/
    DEST_DIR=$OPA_HOME/wrkdir/MODEL/DA__FREQ_1/
 OUTNC_CHECK=$OPA_HOME/wrkdir/float_preproc/OUTNC/
OUTTXT_CHECK=$OPA_HOME/wrkdir/float_preproc/OUTTXT/

export ONLINE_REPO=/gpfs/scratch/userexternal/gbolzon0/CHAIN_V5C/ONLINE_V5C
export MASKFILE=$OPA_HOME/wrkdir/MODEL/meshmask.nc
mkdir -p $OUTNC_CHECK
mkdir -p $OUTTXT_CHECK
mkdir -p $BASEDIR
mkdir -p $DA_DIR/links
mv $DA_DIR/RSTbefore.*0000* $DA_DIR/links

# DADEP=600

for vv in N3n P_l; do

	rm -rf $BASEDIR
	mkdir -p $BASEDIR
	cd $DIR

	if [ $vv == P_l ]; then

	    python var_aggregator.py -l RSTbefore.${DATE}*13:00*P1l.nc -i $DA_DIR -d VarDescriptor_P_lagg.xml -t $TMP_DIR -m $MASKFILE
	    MODEL_AVEDIR=$TMP_DIR
	    DADEP=200

	fi

	if [ $vv == N3n ] ; then
	    MODEL_AVEDIR=$DA_DIR
	    DADEP=600
	fi

	VAR_DESCRIPTOR=$OPA_HOME/wrkdir/float_preproc/VarDescriptor_${vv}.xml
	python float_extractor.py -t ${DATE}  -i $MODEL_AVEDIR -b $BASEDIR -d $VAR_DESCRIPTOR -v $vv

	python preproc.py              -t ${DATE}  -i $MODEL_AVEDIR -b $BASEDIR -m $MASKFILE -v $vv -d $DADEP



done

mv $DA_DIR/links/RSTbefore* $DA_DIR

echo merge arg_mis and copy
echo ---  python merge_arg_mis.py -t ${DATE}
python merge_arg_mis.py -t ${DATE}
mv ${DATE}.arg_mis.dat $DEST_DIR


for vv in N3n P_l; do

mv ${DATE}.${vv}_arg_mis.dat $DEST_DIR

done




