#!/bin/bash

usage() {
echo "SYNOPSYS"
echo "config_OGSTM.sh arch "
}

if [ $# -lt 1 ] ; then
   usage
   exit 1
fi

OPA_ARCH=$1
#OMP_FLAG=$2
export OGSTMDIR=$PWD/..



#----------------- User configuration -----------------
archfile="${OGSTMDIR}/compilers/compiler.inc"
if [ $OPA_ARCH == 'ppc64' ] ; then
   cppdefs="-WF,-Dkey_trahdfcoef1d,-Dkey_trahdfbilap,-Dkey_mpp,-Dkey_mpp_mpi,-Dkey_trc_smolar,-Dkey_trc_hdfbilap,-Dkey_trc_dmp,-Dkey_kef,-Dkey_trc_sed,-Dkey_trc_bfm,-Dkey_INCLUDE_BFM_PELCO2"

   if [ $OPENMP_FLAG ] ; then cppdefs+=",-D__OPENMP"; fi
else
   cppdefs="-Dkey_trahdfcoef1d -Dkey_trahdfbilap -Dkey_mpp -Dkey_mpp_mpi -Dkey_trc_smolar  -Dkey_trc_hdfbilap -Dkey_trc_dmp -Dkey_kef -Dkey_trc_sed -Dkey_trc_bfm -Dkey_INCLUDE_BFM_PELCO2 -DMem_Monitor"
   if [ $OPENMP_FLAG ] ; then cppdefs+=" -D__OPENMP"; fi
fi
exe=${OGSTMDIR}/bin/ogstm.xx
#----------------- User configuration -----------------

if [ "${OGSTMDIR}" = "" ]; then
   echo "Environmental variable OGSTMDIR not defined!"
   exit 1
fi

BLDDIR="${OGSTMDIR}/build/BLD_OGSTM"
MKMF="${OGSTMDIR}/bin/mkmf"
oflags=" -I${OGSTMDIR}/src/ "

if [ ! -d ${BLDDIR} ]; then
  mkdir ${BLDDIR}
fi

cd ${BLDDIR}

ls ${OGSTMDIR}/src/General/get_mem.c       >  OGSTM.lst
ls ${OGSTMDIR}/src/General/get_mem_mod.F90 >> OGSTM.lst
${MKMF} -c "${cppdefs}" -o "${oflags}" -t ${archfile} -p libNotUsed.a OGSTM.lst
mv Makefile MakeLib

ls ${OGSTMDIR}/src/General/*.[Fhc] >  OGSTM.lst
ls ${OGSTMDIR}/src/IO/*.[Fhc]     >> OGSTM.lst
ls ${OGSTMDIR}/src/MPI/*.[Fhc]    >> OGSTM.lst
ls ${OGSTMDIR}/src/PHYS/*.[Fhc]   >> OGSTM.lst
ls ${OGSTMDIR}/src/BIO/*.[Fhc]    >> OGSTM.lst


${MKMF} -c "${cppdefs}" -o "${oflags}" -t ${archfile} -p ${exe} OGSTM.lst
