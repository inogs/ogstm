#! /bin/bash
OGSTM_ARCH=x86_64
OGSTM_OS=LINUX
OGSTM_COMPILER=intel
DEBUG=
DEBUG=.dbg
BFM_PRESET=OGSTM_BENTHIC.B1B
#DEBUG=   # for production
ISFLUXUS=false

export OPENMP_FLAG=     #-fopenmp



usage() {
echo "SYNOPSYS"
echo "Build BFM and ogstm model"
echo "ogstm_bfm_builder.sh [ BFMDIR ] [ OGSTMDIR ]"
echo ""
echo " Dirs have to be expressed as full paths "
echo "EXAMPLE"
echo " ./ogstm_bfm_builder.sh $PWD/bfm $PWD/ogstm "

}

if [ $# -lt 2 ] ; then
   usage
   exit 1
fi

BFMDIR=$1
OGSTMDIR=$2


############### MODULES AND ENVIRONMENT VARIABLES

if [[ $OGSTM_COMPILER == gfortran ]] ; then
   if [ $ISFLUXUS == true ] ; then
       module purge
       module load openmpi-x86_64
   else
       module load autoload openmpi/1.4.4--gnu--4.5.2 netcdf/4.1.3--gnu--4.5.2 #plx
   fi

else
   if [[ $OGSTM_ARCH == ppc64 ]] ; then
       module load profile/advanced
       module load bgq-xl/1.0  netcdf/4.1.3--bgq-xl--1.0
       module load hdf5/1.8.9_ser--bgq-xl--1.0 szip/2.1--bgq-xl--1.0 zlib/1.2.7--bgq-gnu--4.4.6
   fi
   if [[ $OGSTM_ARCH == x86_64 ]] ; then
      if  [ $ISFLUXUS == true ] ; then
         module purge
         module load intel-openmpi
         export NETCDF_INC=/usr/local/intel/include
         export NETCDF_LIB=/usr/local/intel/lib

      else
         module load profile/advanced
         module load intel/cs-xe-2015--binary intelmpi/5.0.1--binary
         module load hdf5/1.8.13_ser--intel--cs-xe-2015--binary netcdf/4.1.3--intel--cs-xe-2015--binary
      fi
   fi
fi

##############################################################



# ----------- BFM library ---------------------
cd $BFMDIR
A=`svn info 2>&1 `  # exit status 1 if bfmv5
if [ $? == 0 ] ; then 
   BFMversion=BFMv2
else
   BFMversion=bfmv5
fi

INC_FILE=${OGSTM_ARCH}.${OGSTM_OS}.${OGSTM_COMPILER}${DEBUG}.inc

if [ $BFMversion == BFMv2 ] ; then
   cd ${BFMDIR}/compilers
   cp $INC_FILE compiler.inc

   ###################################################### just because R1.3 does not have include/
   mkdir -p  ${BFMDIR}/include

   cd ${BFMDIR}/build
   ./config_BFM.sh -a ${OGSTM_ARCH} -c ogstm
   cd BLD_OGSTMBFM
   gmake

else
   # in-place replace the entire ARCH line
   sed -i "s/.*ARCH.*/        ARCH    = '$INC_FILE'  /"  build/configurations/OGS_PELAGIC/configuration
   cd $BFMDIR/build
   ./bfm_configure.sh -gc -o ../lib/libbfm.a -p $BFM_PRESET
   if [ $? -ne 0 ] ; then  echo  ERROR; exit 1 ; fi
fi

export BFM_INC=${BFMDIR}/include
export BFM_LIB=${BFMDIR}/lib




# ------------ OGSTM builder -----------------

cd ${OGSTMDIR}/compilers
INC_FILE=${OGSTM_ARCH}.${OGSTM_OS}.${OGSTM_COMPILER}${DEBUG}.inc
cp $INC_FILE compiler.inc

cd ${OGSTMDIR}/build
./config_OGSTM.sh ${OGSTM_ARCH} 
cd BLD_OGSTM
make -f MakeLib
rm -f get_mem_mod.o
gmake

if [ $? -ne 0 ] ; then  echo  ERROR; exit 1 ; fi

### OGSTM NAMELIST GENERATION (also by Frequency Control )

mkdir -p ${OGSTMDIR}/ready_for_model_namelists/

if [ $BFMversion == bfmv5 ] ; then
   cp ${BFMDIR}/build/tmp/${BFM_PRESET}/namelist.passivetrc ${OGSTMDIR}/bfmv5/
   cd ${OGSTMDIR}/bfmv5/
   python ogstm_namelist_gen.py #generates namelist.passivetrc_new

   cp ${OGSTMDIR}/src/namelists/namelist*      ${OGSTMDIR}/ready_for_model_namelists/ 
   cp namelist.passivetrc_new                  ${OGSTMDIR}/ready_for_model_namelists/namelist.passivetrc #overwriting
   cp ${BFMDIR}/build/tmp/${BFM_PRESET}/*.nml  ${OGSTMDIR}/ready_for_model_namelists/
else
   #V2
   cp ${OGSTMDIR}/src/namelists/namelist*    ${OGSTMDIR}/ready_for_model_namelists/
   cp ${BFMDIR}/src/namelist/*.nml           ${OGSTMDIR}/ready_for_model_namelists/
fi
