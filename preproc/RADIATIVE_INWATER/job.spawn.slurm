#!/bin/bash

#SBATCH --job-name=POST
#SBATCH -N1
#SBATCH --ntasks-per-node=36
#SBATCH --time=6:30:00
#SBATCH --mem=100gb
#SBATCH --account=OGS_dev_1
#SBATCH --partition=gll_usr_prod

cd $SLURM_SUBMIT_DIR

source ~/sequence.sh
module load profile/base
module load autoload
module load hdf5/1.8.18--intel--pe-xe-2018--binary netcdf/4.6.1--intel--pe-xe-2018--binary
module load mpi4py/3.0.0--intelmpi--2018--binary
export PYTHONPATH=/cineca/prod/opt/libraries/numpy/1.15.2/python--2.7.12/lib/python2.7/site-packages/:
/galileo/home/userexternal/eterzic0/bit.sea


source /gpfs/work/OGS18_PRACE_P_0/COPERNICUS/py_env_2.7.12/bin/activate

module load profile/advanced
module load intel/pe-xe-2018--binary
module load intelmpi/2018--binary
module load netcdf/4.6.1--intel--pe-xe-2018--binary
module load netcdff/4.4.4--intel--pe-xe-2018--binary
module load cmake/3.12.0
module load numpy/1.15.2--python--2.7.12

date
ipython launcher_parall.py
date


