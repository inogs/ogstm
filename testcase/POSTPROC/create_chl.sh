#! /bin/bash

cp -f ave.P1l.nc ave.aux1.nc
cp -f ave.P2l.nc ave.aux2.nc
cp -f ave.P3l.nc ave.aux3.nc
cp -f ave.P4l.nc ave.aux4.nc

ncrename -v P1l,chl ave.aux1.nc
ncrename -v P2l,chl ave.aux2.nc
ncrename -v P3l,chl ave.aux3.nc
ncrename -v P4l,chl ave.aux4.nc

ncbo --op_typ=add ave.aux1.nc ave.aux2.nc ave.aux5.nc
ncbo --op_typ=add ave.aux5.nc ave.aux3.nc ave.aux6.nc
ncbo --op_typ=add ave.aux6.nc ave.aux4.nc ave.chl.nc
rm -f ave.aux?.nc
