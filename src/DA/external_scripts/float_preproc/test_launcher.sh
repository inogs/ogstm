#! /bin/bash

module purge
module load profile/base
module load intel/pe-xe-2018--binary intelmpi/2018--binary
module load autoload
module load hdf5/1.8.18--intel--pe-xe-2018--binary netcdf/4.6.1--intel--pe-xe-2018--binary
module load numpy/1.15.2--python--2.7.12 #mpi4py/3.0.0--intelmpi--2018--binary
module load ncview


export     ONLINE_REPO=/gpfs/scratch/userexternal/gbolzon0/CHAIN_V5C/ONLINE_V5C
export        MASKFILE=/gpfs/work/OGS_prod_0/OPA/V6C/devel/wrkdir/analysis/2/MODEL/meshmask.nc
export      OPA_SCRDIR=/gpfs/scratch/userexternal/gbolzon0/CHAIN_V6C/bin
export          DA_DIR=/gpfs/scratch/userexternal/gbolzon0/CHAIN_V6C/TEST/DA
export      DA__FREQ_1=/gpfs/scratch/userexternal/gbolzon0/CHAIN_V6C/TEST/DA__FREQ_1
export      OPA_VENV_1=/gpfs/work/OGS20_PRACE_P/COPERNICUS/py_env_2.7.12/
export      OPA_BITSEA=/gpfs/work/OGS20_PRACE_P/COPERNICUS/bit.sea
export      OPA_BINDIR=/gpfs/scratch/userexternal/gbolzon0/CHAIN_V6C/CODE/ogstm/src/DA/external_scripts/float_preproc
export OPA_POSTPROCDIR=/gpfs/scratch/userexternal/gbolzon0/CHAIN_V6C/bin


./Float_misfit_gen.sh -t 20190919

