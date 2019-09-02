#!/bin/sh
for AFLE in opt cld aot modcld modcldclim modaer modaerclim aotclim
do 
 if [ -f $AFLE.dat ] ; then rm $AFLE.dat ; fi
done
if [ -f Next.month ] ; then rm Next.month ; fi
