#! /bin/bash

# genInputsDatelists.sh 

# launched from MODEL directory, creates datelists files
# from FORCINGS and BC
# Authomatic generations, no user choice expected.

# Author: GB, 12.03.2012



#################### Function definitions #######################

function datelist {
VAR=$1
for FILE in `ls $VAR*.nc` ; do 
   V=${FILE#$VAR} 
   echo ${V%.nc}  
done


}

function datelist2 {
DIR=$1
VAR=$2
for yyyy in $(ls $DIR); do
  [[  ${#yyyy} == 4 ]] || continue
  for mm in $(ls $DIR/$yyyy) ; do
     [[  ${#mm} == 2 ]] || continue
     STR=$DIR/$yyyy/$mm/${VAR}
     for FILE in $(ls ${STR}* 2>/dev/null ) ; do
           V=${FILE#$STR}
           echo ${V%.nc}
    done
  done
done

}

#################### Input parsing ##############################

if [ $# -gt 1 ]; then BCS_FILE=$1; else BCS_FILE="boundaries.nml"; fi


#################### O-O boundary conditions ####################

# Compatibility with old BCS names
declare -A BCS_NAMES
BCS_NAMES=(
    ["riv"]="TIN"
    ["atl"]="ATL"
    ["gi1"]="GIB"
    ["gi2"]="GIB"
    ["gi3"]="GIB"
)

tail -n +2 $BCS_FILE | while read BC_raw; do

    BC=$(echo $BC_raw | sed 's/"//g')

    BC_NAME=$(echo $BC | cut -d ' ' -f 1 | sed 's/,//g'); echo "BC_NAME="$BC_NAME
    BC_LIST=$(echo $BC | cut -d ' ' -f 4 | sed 's/,//g')
    BC_PERIODIC=$(echo $BC | cut -d ' ' -f 5 | sed 's/,//g')

    if [ $BC_PERIODIC == "T" ]; then
	  echo "BC/${BCS_NAMES[$BC_NAME]}_yyyy*"
        ls BC/${BCS_NAMES[$BC_NAME]}_yyyy* | wc -l > $BC_LIST
        ls BC/${BCS_NAMES[$BC_NAME]}_yyyy* | sed "s/^/'/g" | sed "s/$/'/g" >> $BC_LIST
    else
        ls BC/${BCS_NAMES[$BC_NAME]}_[0-9]* | wc -l > $BC_LIST
        ls BC/${BCS_NAMES[$BC_NAME]}_[0-9]* | sed "s/^/'/g" | sed "s/$/'/g" >> $BC_LIST
    fi

done



#################### Previous script ############################

datelist BC/ATM_      > AtmTimes
datelist BC/CO2_      > carbonTimes

datelist2 FORCINGS U   > forcingsTimes
if [ -d OPTICS ] ; then
   datelist2 OPTICS atm.     >  optatmTimes
   datelist OPTICS/climatm.  >  optclimTimes
   datelist OPTICS/aero.     >  optAeroTimes
else
  datelist KEXT/KextF_  > kextTimes
fi

#for I in `ls SATELLITE/ | cut -c 1-8` ; do echo $I-12:00:00 ; done > daTimes_sat
