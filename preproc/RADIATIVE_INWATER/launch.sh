#! /bin/bash

source compila_MARCONI.sh

while read yyyymmdd ;
do

    echo "$yyyymmdd"

    INPUT_OASIM_FILE='/marconi_work/OGS_dev_0/BIOPTIMOD_DAILY/oasimpub/data/'${yyyymmdd}'/rad_0m'${yyyymmdd}'.nc'
    INPUT_ACDOM_FILE='boussole_absorption/aCDOM_boussole_'${yyyymmdd}'.nc'
    INPUT_APHY_FILE='boussole_absorption/aphy_boussole_'${yyyymmdd}'.nc'
    INPUT_ANAP_FILE='boussole_absorption/anap_boussole_'${yyyymmdd}'.nc'
    INPUT_PFT_FILE='pigments/PFT_boussole_'${yyyymmdd}'.nc'

    OUTPUT_FILE='output/edeu_boussole_'${yyyymmdd}'.nc'

    ./compute.xx $INPUT_OASIM_FILE $INPUT_ACDOM_FILE $INPUT_APHY_FILE $INPUT_ANAP_FILE $INPUT_PFT_FILE $OUTPUT_FILE


done < lista_date.txt
