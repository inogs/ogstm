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
export ONLINE_REPO=/gpfs/scratch/userexternal/gbolzon0/CHAIN_V5C/ONLINE_V5C
export MASKFILE=$OPA_HOME/wrkdir/MODEL/meshmask.nc


     CODEDIR=$OPA_HOME/wrkdir/float_preproc/
      DA_DIR=$OPA_HOME/wrkdir/MODEL/DA__FREQ_1/
     BASEDIR=$DA_DIR/PROFILATORE
    DEST_DIR=$DA_DIR


##  V6C version
#. @@(I:OPA_HOME)/bin/opa_profile.inc
# opa_prex "module unload numpy"
# opa_prex "source $OPA_VENV_1/bin/activate"
# PYTHONPATH=${PYTHONPATH}:$OPA_BITSEA
# DA_DIR=$OPA_WRKDIR/MODEL/DA__FREQ_1/
# BASEDIR=$DA_DIR/PROFILATORE
# export MASKFILE=$OPA_WRKDIR/MODEL/meshmask.nc
# export ONLINE_REPO=$I_OPA_HOME/inpdir/ONLINE


mkdir -p $DA_DIR/links
mv $DA_DIR/RSTbefore.*0000* $DA_DIR/links


for vv in N3n P_l; do

	rm -rf $BASEDIR
	mkdir -p $BASEDIR

    VAR_DESCRIPTOR="$OPA_HOME/wrkdir/float_preproc/VarDescriptor_${vv}.xml" # esterni
    MISFIT_FILE="${DA_DIR}/${DATE}.${vv}_arg_mis.dat"

	if [ $vv == P_l ]; then
	    python ${CODEDIR}/var_aggregator.py -l RSTbefore.${DATE}*13:00*P1l.nc -i $DA_DIR -d $VAR_DESCRIPTOR -t $DA_DIR -m $MASKFILE
	    DADEP=200
	fi

	if [ $vv == N3n ] ; then
	    DADEP=600
	fi


	python ${CODEDIR}/float_extractor.py -t ${DATE}  -i $DA_DIR -b $BASEDIR -d $VAR_DESCRIPTOR -v $vv
	python ${CODEDIR}/preproc.py  -t ${DATE}  -i $DA_DIR -b $BASEDIR -m $MASKFILE -v $vv -d $DADEP --misfit $MISFIT_FILE -o $OUTDIR

done

python ${CODEDIR}/merge_arg_mis.py -n ${DA_DIR}/${DATE}.N3n_arg_mis.dat -c ${DA_DIR}/${DATE}.P_l_arg_mis.dat -o ${DA_DIR}/${DATE}.arg_mis.dat


mv $DA_DIR/links/RSTbefore* $DA_DIR




