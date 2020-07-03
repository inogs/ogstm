#! /bin/bash

source ~/sequence.sh

while read yyyymmdd ;
do

echo "$yyyymmdd"

python create_ave.py -d ${yyyymmdd}

done < lista_date.txt