#! /bin/bash

HERE=$PWD
ORIG_DIR=/pico/scratch/userexternal/gbolzon0/Carbonatic-01/FORCINGS/PHYS/UNZIPPED
cd $ORIG_DIR
DATE=2015
#DATE=20150[56]
ls -1 $DATE*U.nc > $HERE/nomefile_U
ls -1 $DATE*V.nc > $HERE/nomefile_V
ls -1 $DATE*W.nc > $HERE/nomefile_W
ls -1 $DATE*T.nc > $HERE/nomefile_T

