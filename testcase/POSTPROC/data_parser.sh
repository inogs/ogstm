#! /bin/bash

module load autoload nco

TEST=$1
MODDIR="../${TEST}/wrkdir/MODEL"

ncrcat -O -h ${MODDIR}/AVE_FREQ_1/ave.????????-??:??:??.nc aux.nc
ncks -O --mk_rec_dmn time aux.nc ave.${TEST}.nc
rm aux.nc

