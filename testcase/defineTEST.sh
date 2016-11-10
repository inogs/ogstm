#! /bin/bash

for area in  ALBORAN DYFAMED LEVANTINE NA ALBORAN_ATM DYFAMED_ATM LEVANTINE_ATM DYFAMED_DIA; 
do
   echo $area
 cp -f TEST01-REF/wrkdir/MODEL/1.aveTimes TEST_${area}/wrkdir/MODEL/
 cp -f TEST01-REF/wrkdir/MODEL/Start_End_Times TEST_${area}/wrkdir/MODEL/
 cp -f TEST01-REF/wrkdir/MODEL/namelist.init TEST_${area}/wrkdir/MODEL/
 cp -f TEST01-REF/wrkdir/MODEL/job.pbs TEST_${area}/wrkdir/MODEL/
 cp -f FLUXES/Fluxes.nc TEST_${area}/wrkdir/MODEL/
done
