#! /bin/bash
  
#SBATCH --job-name=ESIWACE_TEST
#SBATCH --account=OGS21_PRACE_P_0
#SBATCH --partition=m100_usr_prod
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=32
#SBATCH --gres=gpu:4
#SBATCH --constraint=gpureport
#SBATCH --exclusive


source $CINECA_SCRATCH/esiwace_m100/ModelBuild_mpi/ogstm/compilers/machine_modules/m100.hpc-sdk
export RANKS_PER_NODE=$SLURM_TASKS_PER_NODE
export SRC=$(realpath $1)
export DST=$TMPDIR/ESIWACE
export WRK=${DST}/$(basename ${SRC})


mkdir -p ${DST}
cp -r "${SRC}" "${DST}"
cd "${WRK}" || exit 

mpirun -gpu -np $RANKS_PER_NODE ./ogstm.xx

