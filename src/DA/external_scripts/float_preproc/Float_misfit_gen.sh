#! /bin/bash

usage() {
echo "Generates misfit files for argo floats "
echo "SYNOPSYS"
echo "Float_misfit_gen.sh [ -t date] "
echo "EXAMPLE"
echo 'Float_misfit_gen.sh -t 20150106 '
echo ""
}


if [ $# -lt 2 ] ; then
  usage
  exit 1
fi

for I in 1; do
   case $1 in
      "-t" ) DATE=$2;;
        *  ) echo "Unrecognized option $1." ; usage;  exit 1;;
   esac
   shift 2
done


module purge
module load profile/base
module load autoload
module load netcdf/4.6.1--intel--pe-xe-2018--binary
module load numpy/1.15.2--python--2.7.12 
module load intelmpi/2018--binary

source /gpfs/work/OGS20_PRACE_P/COPERNICUS/py_env_2.7.12/bin/activate
export PYTHONPATH=$PYTHONPATH:/gpfs/work/OGS20_PRACE_P/COPERNICUS/bit.sea

export OPA_HOME=/gpfs/scratch/userexternal/ateruzzi/MULTIVARIATE_24/TEST_04/

     CODEDIR=$OPA_HOME/wrkdir/float_preproc/
      DA_DIR=$OPA_HOME/wrkdir/MODEL/DA__FREQ_1/
     BASEDIR=$DA_DIR/PROFILATORE
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



for vv in N3n P_l; do

	rm -rf $BASEDIR
	mkdir -p $BASEDIR


	if [ $vv == P_l ]; then

	    python ${CODEDIR}/var_aggregator.py -l RSTbefore.${DATE}*13:00*P1l.nc -i $DA_DIR -d ${CODEDIR}/VarDescriptor_P_l.xml -t $DA_DIR -m $MASKFILE
	    DADEP=200

	fi

	if [ $vv == N3n ] ; then
	    DADEP=600
	fi

	VAR_DESCRIPTOR=$OPA_HOME/wrkdir/float_preproc/VarDescriptor_${vv}.xml
	python ${CODEDIR}/float_extractor.py -t ${DATE}  -i $DA_DIR -b $BASEDIR -d $VAR_DESCRIPTOR -v $vv

	python ${CODEDIR}/preproc.py   -t ${DATE}  -i $DA_DIR -b $BASEDIR -m $MASKFILE -v $vv -d $DADEP



done

mv $DA_DIR/links/RSTbefore* $DA_DIR


python ${CODEDIR}/merge_arg_mis.py -t ${DATE}

mv ${DATE}.arg_mis.dat $DEST_DIR


for vv in N3n P_l; do
    mv ${DATE}.${vv}_arg_mis.dat $DEST_DIR
done




