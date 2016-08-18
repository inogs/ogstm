#! /bin/bash
for area in  ALBORAN DYFAMED LEVANTINE NA; 
do
   echo $area
   mkdir -p TEST_${area}/wrkdir/MODEL
   rsync -avz -e ssh plazzari@login.pico.cineca.it:/pico/scratch/userexternal/plazzari/TILMAN/ogstm/testcase/TEST_${area}/wrkdir/MODEL/meshmask.nc TEST_${area}/wrkdir/MODEL/meshmask.nc
done

rsync -avz -e ssh plazzari@login.pico.cineca.it:/pico/scratch/userexternal/plazzari/TILMAN/ogstm/testcase/POSTPROC .
rsync -avz -e ssh plazzari@login.pico.cineca.it:/pico/scratch/userexternal/plazzari/TILMAN/ogstm/testcase/TEST_LIST.dat .
