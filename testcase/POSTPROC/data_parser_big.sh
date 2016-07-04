#! /bin/bash

module load autoload nco

TEST='TEST02'
MODDIR="../${TEST}/wrkdir/MODEL"
VAR=$1

ncrcat -O -h ${MODDIR}/AVE_FREQ_1/ave.????????-??:??:??.nc aux.nc
ncks -O --mk_rec_dmn time aux.nc ave.${VAR}.nc
rm aux.nc


