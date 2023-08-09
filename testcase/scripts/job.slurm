#! /bin/bash
  
#SBATCH --job-name=OGSTM-BFM_GPU_TEST
#SBATCH --account=OGS23_PRACE_IT_0
#SBATCH --partition=boost_usr_prod
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:4
#SBATCH --exclusive

source "~/ModelBuild/ogstm/compilers/machine_modules/leonardo.nvhpc"
export RANKS_PER_NODE=$SLURM_TASKS_PER_NODE
export SRC=$(realpath $1)
export DST=$TMPDIR/OGSTM-BFM_GPU_TEST
export WRK=${DST}/$(basename ${SRC})

mkdir -p ${DST}
cp -r "${SRC}" "${DST}"
cd "${WRK}" || exit 

mpirun -gpu -np $RANKS_PER_NODE ./ogstm.xx

