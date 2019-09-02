#! /bin/bash
#
#SBATCH -N 1 -n 1              # 1 tasks on 1 node
#SBATCH --time=4:00:00         # time limits: 20 min
#SBATCH --error myJob.err      # std-error file
#SBATCH --output myJob.out     # std-output file
#SBATCH --account=OGS18_PRACE_P  # account number
#SBATCH --partition=bdw_all_serial  # partition to be used

#

cd $SLURM_SUBMIT_DIR
date

source $HOME/sequence.sh

for year in 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016 2017; do

    python mars_request_daily_1YEAR.py $year

done

date
