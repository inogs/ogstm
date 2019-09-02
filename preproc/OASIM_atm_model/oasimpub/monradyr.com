#!/bin/bash  
IRRDIR=/gpfs/scratch/userexternal/eterzic0/BIOPTIMOD_HF/DATA
yyyymmdd=`cat bcs/yr.dat`
YR=$( echo $yyyymmdd | cut -c 1-4)
MON=$( echo $yyyymmdd | cut -c 5-6)
DAY=$( echo $yyyymmdd | cut -c 7-8)
JDAY=$(echo $(date --date=$yyyymmdd +%j))
echo ' yyyymmdd ' $yyyymmdd
echo ' YR ' $YR
echo ' MON ' $MON
echo ' DAY ' $DAY
echo 'JDAY ' $JDAY
IRRDIRYR=$IRRDIR/$yyyymmdd
if [ ! -d $IRRDIRYR ] ; then mkdir $IRRDIRYR ; fi
cd src
rm *.o
#make monrad
./compile.sh
cp monrad ../
cd ../
#  
# Set up input data files
for AFLE in opt cld aot atmo25b abw25b acbc25b clddays slingo
do
 if [ -f bcs/$AFLE.dat.gz ]
 then
  echo 'Uncompressing ' $AFLE.dat.gz
  gunzip bcs/$AFLE.dat.gz
 fi
done
cp bcs/opt.dat opt.dat
cp bcs/cld.dat cld.dat

ln -fs /gpfs/scratch/userexternal/plazzari/BIOPTIMOD_HF/tools/CREATE_INPUT_ECMWF/output/opt${yyyymmdd}_ECMWF.nc opt.nc
ln -fs /gpfs/scratch/userexternal/plazzari/BIOPTIMOD_HF/tools/CREATE_INPUT_ECMWF/output/clouds${yyyymmdd}_ECMWF.nc clouds.nc

if [ -f modcld.dat ] ; then rm modcld.dat ; fi
if [ -f atmdata/modcld$YR.dat.gz ] 
then 
 cp atmdata/modcld$YR.dat.gz .
 gunzip modcld$YR.dat.gz
 cp modcld$YR.dat modcld.dat
 rm modcld$YR.dat
fi
cp atmdata/modcld0000.dat.gz modcldclim.dat.gz
gunzip modcldclim.dat.gz
cp atmdata/aot0000.dat.gz aotclim.dat.gz
gunzip aotclim.dat.gz
if [ -f bcs/aot.dat ] ; then cp bcs/aot.dat aot.dat ; fi
#
#  Run
echo $JDAY > day.dat
echo $MON > month.dat
echo $YR  > year.dat
if [ ! -f $IRRDIRYR/swr$yyyymmdd.dat.gz ]
 then
  if [ -f modaer.dat ] ; then rm modaer.dat ; fi
  if [ -f modaerclim.dat ] ; then rm modaerclim.dat ; fi
  if [ -f atmdata/modaer$YR$MON.dat.gz ]
  then
   cp -p atmdata/modaer$YR$MON.dat.gz modaer.dat.gz
   gunzip modaer.dat.gz
  else
   cp -p atmdata/modaer$MON.dat.gz modaerclim.dat.gz
   gunzip modaerclim.dat.gz
  fi
  ./monrad > out
  mv out $IRRDIRYR/out$yyyymmdd
  if [ -f modaer.dat ] ; then rm modaer.dat ; fi
  if [ -f modaerclim.dat ] ; then rm modaerclim.dat ; fi

  mv rad rad$yyyymmdd.dat
  gzip rad$yyyymmdd.dat
  mv rad$yyyymmdd.dat.gz $IRRDIRYR
  mv rad_0p.nc rad_0p$yyyymmdd.nc
  mv rad_0m.nc rad_0m$yyyymmdd.nc
  gzip rad_0p$yyyymmdd.nc
  gzip rad_0m$yyyymmdd.nc
  mv rad_0p$yyyymmdd.nc.gz $IRRDIRYR
  mv rad_0m$yyyymmdd.nc.gz $IRRDIRYR

  mv avgirr.dat swr$yyyymmdd.dat
  gzip swr$yyyymmdd.dat
  mv swr$yyyymmdd.dat.gz $IRRDIRYR
  mv swr.nc swr$yyyymmdd.nc
  gzip swr$yyyymmdd.nc
  mv swr$yyyymmdd.nc.gz $IRRDIRYR

  mv avgparq.dat parq$yyyymmdd.dat
  gzip parq$yyyymmdd.dat
  mv parq$yyyymmdd.dat.gz $IRRDIRYR
  mv parq.nc parq$yyyymmdd.nc
  gzip parq$yyyymmdd.nc
  mv parq$yyyymmdd.nc.gz $IRRDIRYR

  mv eds.dat eds$yyyymmdd.dat
  gzip eds$yyyymmdd.dat
  mv eds$yyyymmdd.dat.gz $IRRDIRYR
  mv eds.nc eds$yyyymmdd.nc
  gzip eds$yyyymmdd.nc
  mv eds$yyyymmdd.nc.gz $IRRDIRYR

fi
#
# Remove  input data files
for AFLE in modcld opt cld aot aotclim modcldclim
do
 if [ -f bcs/$AFLE.dat ]
 then
  rm bcs/$AFLE.dat
 fi
 if [ -f $AFLE.dat ]
 then
  rm $AFLE.dat
 fi
done
