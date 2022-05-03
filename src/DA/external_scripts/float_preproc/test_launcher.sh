#! /bin/bash

module purge
module load autoload
module load autoload
module load intel/oneapi-2021--binary
module load intelmpi/oneapi-2021--binary
#module load numpy/1.15.2--python--2.7.12 #mpi4py/3.0.0--intelmpi--2018--binary


BASE=/g100_scratch/userexternal/camadio0
export     ONLINE_REPO=/g100_scratch/userexternal/camadio0/SUPERFLOAT_2012_2021_V8c/
export        MASKFILE=$BASE/Multivariate_Assimilation_2017_2018/wrkdir/MODEL/meshmask.nc
export      OPA_SCRDIR=$BASE/Multivariate_Assimilation_2017_2018/wrkdir/POSTPROC/bin
export          DA_DIR=$BASE/Multivariate_Assimilation_2017_2018/wrkdir/DA
export      DA__FREQ_1=$BASE/Multivariate_Assimilation_2017_2018/wrkdir/MODEL/DA__FREQ_1
export      OPA_VENV_1=/g100_work/OGS20_PRACE_P_2/COPERNICUS/py_env_3.6.8/
export      OPA_BITSEA=/g100/home/userexternal/camadio0/bit.sea_py3/
export      OPA_BINDIR=$BASE/Multivariate_Assimilation_2017_2018/wrkdir/float_preproc
export OPA_POSTPROCDIR=$BASE/Multivariate_Assimilation_2017_2018/wrkdir/POSTPROC/bin

# edit Float_misfit_gen.sh
# rename opa_prex_ -->  my_prex
# rename opa_profile --> profile
# comment module unload numpy, is useful only for chain 

mv end.txt start.txt
./Float_misfit_gen.sh #-t 20170101


