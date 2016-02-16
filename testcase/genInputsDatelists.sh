#! /bin/bash

# genInputsDatelists.sh 

# launched from MODEL directory, creates datelists files
# from FORCINGS and BC
# Authomatic generations, no user choice expected.

# Author: GB, 12.03.2012

function datelist {
VAR=$1
for FILE in `ls $VAR*.nc` ; do 
   V=${FILE#$VAR} 
   echo ${V%.nc}  
done


}

DIR=$1	

datelist ${DIR}/BC/TIN_          > ${DIR}/RiversTimes
datelist ${DIR}/BC/GIB_          > ${DIR}/GibTimes
datelist ${DIR}/BC/ATM_          > ${DIR}/AtmTimes
datelist ${DIR}/BC/CO2_          > ${DIR}/carbonTimes

datelist ${DIR}/FORCINGS/U       > ${DIR}/forcingsTimes
datelist ${DIR}/FORCINGS/KextF_  > ${DIR}/kextTimes
