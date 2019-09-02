#! /bin/bash
#
#SBATCH -N 1 -n 1              # 1 tasks on 1 node
#SBATCH --time=24:00:00         # time limits: 20 min
#SBATCH --error myJob.err      # std-error file
#SBATCH --output myJob.out     # std-output file
#SBATCH --account=OGS18_PRACE_P_0  # account number
#SBATCH --partition=gll_usr_prod  # partition to be used


#

cd $SLURM_SUBMIT_DIR

date

source $HOME/sequence.sh

./launcher.sh

date
