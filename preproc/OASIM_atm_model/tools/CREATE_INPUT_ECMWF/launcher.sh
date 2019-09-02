#! /bin/bash

source ~/sequence.sh

while read yyyymmdd ;
do

echo "$yyyymmdd"

python  create_opt_nc.py    ${yyyymmdd}
python  create_cloud_nc.py  ${yyyymmdd}

done < lista_date.txt
