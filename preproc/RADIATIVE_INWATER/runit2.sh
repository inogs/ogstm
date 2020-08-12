#!/bin/bash

##SBATCH --job-name=1D-RTM
#SBATCH -N1
#SBATCH --ntasks-per-node=36
#SBATCH --time=0:30:00 
##SBATCH --mem=100gb    
#SBATCH --account=OGS_dev_1
#SBATCH --partition=gll_usr_prod

cd $SLURM_SUBMIT_DIR

source env.sh

date
srun -n $SLURM_NPROCS python launcher.py
date

mv MATCHUP ../
mv KD ../
mv RRS ../