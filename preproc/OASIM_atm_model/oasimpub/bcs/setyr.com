#!/bin/sh
if [ $# -lt 1 ]
then
 echo 'Enter year to setup'
 exit 1
fi
YR=$1
if [ -f opt.dat.gz ] ; then rm opt.dat.gz ; fi
if [ -f cld1*.dat ] ; then rm cld1*.dat* ; fi
if [ -f cld2*.dat ] ; then rm cld2*.dat* ; fi
if [ -f cld.dat.gz ] ; then rm cld.dat.gz ; fi
if [ -f aot.dat.gz ] ; then rm aot.dat.gz ; fi
if [ -f yr.dat ] ; then rm yr.dat ; fi
cp -p  ../atmdata/opt$YR.dat.gz .
gunzip opt$YR.dat.gz
cp -p  opt$YR.dat opt.dat
rm opt$YR.dat
gzip opt.dat
cp -p  ../atmdata/cld$YR.dat.gz .
gunzip cld$YR.dat.gz
cp -p  cld$YR.dat cld.dat
rm cld$YR.dat
gzip cld.dat
if [ -f ../atmdata/aot$YR.dat.gz ]
then
 cp -p  ../atmdata/aot$YR.dat.gz .
 gunzip aot$YR.dat.gz
 cp -p aot$YR.dat aot.dat
 rm aot$YR.dat
 gzip aot.dat
fi
echo $YR >yr.dat
