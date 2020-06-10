#!/bin/bash

#SBATCH --job-name=POST
#SBATCH -N1
#SBATCH --ntasks-per-node=36
#SBATCH --time=6:30:00
#SBATCH --mem=100gb
#SBATCH --account=OGS_dev_1
#SBATCH --partition=gll_usr_prod

cd $SLURM_SUBMIT_DIR

source /galileo/home/userexternal/eterzic0/sequence_eterzic.sh

date
ipython launcher_parall.py
date


