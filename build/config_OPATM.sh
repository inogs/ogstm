#!/bin/sh
export OPATMDIR=$PWD/..
#----------------- User configuration -----------------
archfile="${OPATMDIR}/compilers/compiler.inc"
cppdefs="-Dkey_off_tra -Dkey_trahdfcoef1d -Dkey_trahdfbilap -Dkey_mpp -Dkey_mpp_mpi -Dkey_passivetrc -Dkey_trc_smolar  -Dkey_trc_hdfbilap -Dkey_trc_dmp -Dkey_kef -Dkey_trc_sed -Dkey_trc_bfm -Dkey_INCLUDE_BFM_PELCO2 -D__OPENMP"
exe=${OPATMDIR}/bin/opa.xx
#----------------- User configuration -----------------

if [ "${OPATMDIR}" = "" ]; then
   echo "Environmental variable OPATMDIR not defined!"
   exit 1
fi

BLDDIR="${OPATMDIR}/build/BLD_OPATM"
MKMF="${OPATMDIR}/bin/mkmf"
oflags=" -I${OPATMDIR}/src/ "

if [ ! -d ${BLDDIR} ]; then
  mkdir ${BLDDIR}
fi

cd ${BLDDIR}


ls ${OPATMDIR}/src/General/*.[Fhc] > OPATM.lst
ls ${OPATMDIR}/src/IO/*.[Fhc]     >> OPATM.lst
ls ${OPATMDIR}/src/MPI/*.[Fhc]    >> OPATM.lst
ls ${OPATMDIR}/src/PHYS/*.[Fhc]   >> OPATM.lst
ls ${OPATMDIR}/src/BIO/*.[Fhc]    >> OPATM.lst


${MKMF} -c "${cppdefs}" -o "${oflags}" -t ${archfile} -p ${exe} OPATM.lst
