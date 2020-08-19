#! /bin/bash

source ~/sequence.sh

while read yyyymmdd ;
do

echo "$yyyymmdd"

python create_ave.py ${yyyymmdd}

done < datelist.txt