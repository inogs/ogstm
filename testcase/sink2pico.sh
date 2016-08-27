#! /bin/bash
for area in  ALBORAN_ATM DYFAMED_ATM LEVANTINE_ATM; 
do
   echo $area
   mkdir -p TEST_${area}/wrkdir/MODEL
   rsync -avz -e ssh plazzari@login.pico.cineca.it:/pico/scratch/userexternal/plazzari/TILMAN/ogstm/testcase/TEST_${area}/wrkdir/MODEL/meshmask.nc TEST_${area}/wrkdir/MODEL/meshmask.nc
   rsync -avz -e ssh plazzari@login.pico.cineca.it:/pico/scratch/userexternal/plazzari/TILMAN/ogstm/testcase/TEST_${area}/wrkdir/MODEL/Pelagic_Ecology.nml TEST_${area}/wrkdir/MODEL/Pelagic_Ecology.nml
done

rsync -avz -e ssh plazzari@login.pico.cineca.it:/pico/scratch/userexternal/plazzari/TILMAN/ogstm/testcase/POSTPROC .
rsync -avz -e ssh plazzari@login.pico.cineca.it:/pico/scratch/userexternal/plazzari/TILMAN/ogstm/testcase/POSTPROC_NOATM  .
rsync -avz -e ssh plazzari@login.pico.cineca.it:/pico/scratch/userexternal/plazzari/TILMAN/ogstm/testcase/TEST_LIST.dat .
