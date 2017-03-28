#! /bin/bash

for test in  $(ls -d RUN-??????????); 
do
   echo $test
cp -f 1.aveTimes ${test}/wrkdir/MODEL/
cp -f job.pbs    ${test}/wrkdir/MODEL/
cp -f Fluxes.nc  ${test}/wrkdir/MODEL/
cp -f ${test}/wrkdir/MODEL/Start_End_Times ${test}/wrkdir/MODEL/kextTimes
cd ${test}/wrkdir/MODEL
qsub -N ${test}  job.pbs 
cd ../../../
done
