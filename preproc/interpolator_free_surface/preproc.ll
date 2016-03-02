#!/bin/bash
#
# @ job_name    = FORCING_GEN
#
# @ output      = preproc.$(jobid).out
# @ error       = preproc.$(jobid).err
# @ wall_clock_limit = 6:00:00
# @ class      = serial
# @ job_type   = serial
# @ account_no = IscrA_MUMEBIES
# @ queue




cd /gpfs/scratch/userexternal/ateruzzi/FORCING_RA/INTERP
date

./crea_nomefile.sh
./interp

date

