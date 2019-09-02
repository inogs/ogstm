#!/bin/sh  
IRRDIR=data360180
cd src
rm *.f
rm *.o
make monrad
cp monrad ../
cd ../
#  
# Set up input data files
for afle in cld opt atmo25b abw25b acbc25b clddays \
slingo
do
 if [ -f bcs/$afle.dat.Z ]
 then
  echo 'Uncompressing ' $afle.dat.Z
  uncompress bcs/$afle.dat.Z
 fi
done
cp bcs/cld.dat cld.dat
cp bcs/opt.dat opt.dat
#
#  Run
for MON in 01 02 03 04 05 06 07 08 09 10 11 12
do
 echo $MON >Next.month
 if [ ! -f $IRRDIR/rad$MON.dat.Z ]
 then
  monrad > out
  mv out out$MON
  mv rad $IRRDIR/rad$MON.dat
  compress $IRRDIR/rad$MON.dat
 fi
done
#
#
rm opt.dat
rm cld.dat
rm Next.month
# Compress input data files
for afle in monclimcld monclimopt atmo25b abw25b acbc25b clddays \
slingo
do
 if [ -f bcs/$afle.dat ]
 then
  echo 'Compressing ' $afle.dat
  compress bcs/$afle.dat
 fi
done
