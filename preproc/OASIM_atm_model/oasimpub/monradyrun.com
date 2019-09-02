#!/bin/bash

module purge
module load profile/advanced
module load autoload
module load intel/pe-xe-2018--binary
module load intelmpi/2018--binary
module load netcdf/4.6.1--intel--pe-xe-2018--binary
module load netcdff/4.4.4--intel--pe-xe-2018--binary
module load cmake/3.12.0
module load petsc/3.10.2--intelmpi--2018--binary
module load pnetcdf/1.10.0--intelmpi--2018--binary


WRKDIR=/galileo/home/userexternal/eterzic0/BIOPTIMOD_HF/oasimpub

while read yyyymmdd; do
  echo "$yyyymmdd"

 cd bcs

 ./setyr.com $yyyymmdd

  cd $WRKDIR

 ./monradyr.com >monradyr.doc

done <lista_date.txt
