      parameter(iatm=360,jatm=180)
      parameter(icd=360,jcd=180)
      common /batm/ slp(iatm,jatm),wsm(iatm,jatm),rh(iatm,jatm),
     *oz(iatm,jatm),wv(iatm,jatm),am,Vi
      common /bcld/ ccov(icd,jcd),cldtc(icd,jcd),rlwp(icd,jcd),
     *cdre(icd,jcd)
      common /bcld2/ ndy(icd,jcd),icld(31,31)
      common /baer/ taua(icd,jcd,nlt),asymp(icd,jcd,nlt),
     *ssalb(icd,jcd,nlt)
      real Eda,Esa
      common /beda/ Eda(im,jm,nlt),Esa(im,jm,nlt)
