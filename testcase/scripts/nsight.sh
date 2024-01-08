#! /usr/bin/env bash
#SBATCH --job-name=OGSTM-BFM_NSIGHT_PROFILE
#SBATCH --account=OGS23_PRACE_IT_0
#SBATCH --partition=boost_usr_prod
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:4
#SBATCH --exclusive

module load cuda

ROOT=$CINECA_SCRATCH/ModelBuild

#cd "${ROOT}/ogstm/testcase/$1" || exit

source "${ROOT}/ogstm/compilers/machine_modules/leonardo.nvhpc"

SRC=$(realpath $1)
DST=$TMPDIR/OGSTM-BFM_WORKDIR
WRK=${DST}/$(basename "${SRC}")

mkdir -p "${DST}"
cp -r "${SRC}" "${DST}"
cd "${WRK}" || exit 

RANKS_PER_NODE=$SLURM_TASKS_PER_NODE

rm -rf /tmp/nvidia
ln -s $TMPDIR /tmp/nvidia
nsys profile -o "$PWD/$(date -Iminute)" -f true --stats=true --trace=openacc,mpi mpirun -np $RANKS_PER_NODE ./ogstm.xx
rm -rf /tmp/nvidia
