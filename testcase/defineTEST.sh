#! /bin/bash
#TEST01      TEST01-REF           TEST_ALBORAN  TEST_LEVANTINE  TEST_LIST.dat_template
#TEST01_bad  TEST01_WITH_ATM_RIV  TEST_DYFAMED  TEST_LIST.dat   TEST_NA
for area in  ALBORAN DYFAMED LEVANTINE NA; 
do
   echo $area
 cp -f TEST01-REF/wrkdir/MODEL/1.aveTimes TEST_${area}/wrkdir/MODEL/
 cp -f TEST01-REF/wrkdir/MODEL/Start_End_Times TEST_${area}/wrkdir/MODEL/
 cp -f TEST01-REF/wrkdir/MODEL/namelist.init TEST_${area}/wrkdir/MODEL/
 cp -f TEST01-REF/wrkdir/MODEL/namelist.passivetrc TEST_${area}/wrkdir/MODEL/
 cp -f TEST01-REF/wrkdir/MODEL/job.pbs TEST_${area}/wrkdir/MODEL/
done
