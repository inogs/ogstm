#! /bin/ksh

MPIRUN="mpirun -np 1"

usage() {
echo "Generates misfit files for argo floats "
echo "SYNOPSYS"
echo "Float_misfit_gen.sh"
echo ""
echo "These variables exported by caller:"
echo "ONLINE_REPO"
echo "MASKFILE"
echo "OPA_SCRDIR"
echo "DA_DIR"
echo "DA__FREQ_1"
echo "OPA_VENV_1"
echo "OPA_BITSEA"
echo "OPA_BINDIR"
echo "OPA_POSTPROCDIR"
echo ""
}

echo "Float_misfit_gen.sh start "

source ${OPA_SCRDIR}/opa_profile.inc
#opa_prex "module unload numpy"
opa_prex "source $OPA_VENV_1/bin/activate"
export PYTHONPATH=$OPA_BITSEA:${PYTHONPATH}

while [ True ] ; do
  while [ True ] ; do
    echo "waiting for start.txt"
    # wait of start
    if [ -f start.txt  ] ; then
      break
    fi
    sleep 8
  done


  DATE=`cat start.txt`

  echo "done"


BASEDIR=$DA_DIR/PROFILATORE

mkdir -p ${DA__FREQ_1}/links
mv ${DA__FREQ_1}/RSTbefore.*0000* ${DA__FREQ_1}/links


for vv in N3n P_l O2o; do

	rm -rf $BASEDIR
	mkdir  $BASEDIR

    VAR_DESCRIPTOR="${DA_DIR}/VarDescriptor_${vv}.xml"
       MISFIT_FILE="${DA_DIR}/${DATE}.${vv}_arg_mis.dat"

	if [ $vv == P_l ]; then
           opa_prex_or_die "$MPIRUN python ${OPA_POSTPROCDIR}/var_aggregator.py -l RSTbefore.${DATE}*13:00*P1l.nc -i ${DA__FREQ_1} -d $VAR_DESCRIPTOR -t ${DA__FREQ_1} -m $MASKFILE -s "
            DADEP=200
        fi

        if [ $vv == N3n ] || [ $vv == O2o ]; then
            DADEP=600
        fi


        opa_prex_or_die "$MPIRUN python ${OPA_BINDIR}/float_extractor.py -t ${DATE}  -i ${DA__FREQ_1} -b $BASEDIR -d $VAR_DESCRIPTOR -v $vv "
        opa_prex_or_die "$MPIRUN python ${OPA_BINDIR}/preproc.py  -t ${DATE}  -i ${DA__FREQ_1} -b $BASEDIR -m $MASKFILE -v $vv -d $DADEP --misfit $MISFIT_FILE -o ${DA_DIR} "

done

opa_prex_or_die "python ${OPA_BINDIR}/merge_arg_mis.py -n ${DA_DIR}/${DATE}.N3n_arg_mis.dat -c ${DA_DIR}/${DATE}.P_l_arg_mis.dat -x ${DA_DIR}/${DATE}.O2o_arg_mis.dat -o ${DA__FREQ_1}/${DATE}.arg_mis.dat "


mv ${DA__FREQ_1}/links/RSTbefore* ${DA__FREQ_1}

echo "#############################################"
echo "########### Float_misfit_gen.sh finished ####"
echo "#############################################"

echo "Float_misfit_gen.sh --> dump end.txt"
mv start.txt end.txt

done
