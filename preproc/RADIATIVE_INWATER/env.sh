#!/bin/bash

source ~/sequence.sh

module load profile/base
module load autoload
module load hdf5/1.8.18--intel--pe-xe-2018--binary netcdf/4.6.1--intel--pe-xe-2018--binary
module load mpi4py/3.0.0--intelmpi--2018--binary
module load profile/advanced
module load intel/pe-xe-2018--binary
module load intelmpi/2018--binary
module load netcdf/4.6.1--intel--pe-xe-2018--binary
module load netcdff/4.4.4--intel--pe-xe-2018--binary
module load cmake/3.12.0
module load numpy/1.15.2--python--2.7.12

export PYTHONPATH=/galileo/home/userexternal/eterzic0/BIT_SEA_mod/bit.sea/:$PYTHONPATH

export ONLINE_REPO=/gpfs/scratch/userexternal/eterzic0/FLOAT_BIO/
export STATIC_REPO=/gpfs/scratch/userexternal/eterzic0/Float_OPT_2019/

