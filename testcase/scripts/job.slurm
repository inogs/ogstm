#! /bin/bash
  
#SBATCH --job-name=OGSTM-BFM_GPU_TEST
#SBATCH --account=OGS23_PRACE_IT_0
#SBATCH --partition=boost_usr_prod
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:4
#SBATCH --exclusive

export BUILD_PATH=$(realpath $1)
export SRC=$(realpath $2)
export DST=$(realpath $3)

source "${BUILD_PATH}/machine.modules"
export RANKS_PER_NODE=$SLURM_TASKS_PER_NODE
export WRK="${TMPDIR}/$(basename ${SRC})"

mkdir -p "${WRK}"
cp -r "${SRC}" "${TMPDIR}"
cd "${WRK}" || exit 

ln -sf "${BUILD_PATH}/ogstm.xx"
mpirun -np $RANKS_PER_NODE ./ogstm.xx
cp -r "${WRK}" "${DST}"

