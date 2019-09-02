      subroutine rdatopt(imo,ndmo)
c
c  Reads atmospheric optical data.  
c
#include "define.h"
#include "rad.h"
#include "comlte.h"
#include "comdata.h"
      real ccovm(icd,jcd),cldtcm(icd,jcd)
      real rlwpm(icd,jcd),cdrem(icd,jcd)
      real ccovi(icd,jcd),cldtci(icd,jcd)
      real rlwpi(icd,jcd)
      real slporg(iatm,jatm),wsmorg(iatm,jatm),rhorg(iatm,jatm)
      real ozorg(iatm,jatm),wvorg(iatm,jatm)
      real atmtmp(iatm,jatm)
      real avhrrtau(iatm,jatm),avhrrtauclim(iatm,jatm)
      real facaer(nlt)
      data flag /-0.9E10/
      data flagam /-0.01/
c
      write(6,*)'RDATOPT: imo = ',imo
c  File name: MODIS clouds -- if present use them for re
      ierrmodcld = 0
      imod0 = 0
c  First seek out transient data

      call readnc_global_float_2d('clouds.nc','ccov', ccovm)
      call readnc_global_float_2d('clouds.nc','cldtc', cldtcm)
      call readnc_global_float_2d('clouds.nc','rlwp', rlwpm)
      call readnc_global_float_2d('clouds.nc','cdre', cdrem)

c  Check for valid data
      do j = 1,jcd
       do i = 1,icd
        if (ccovm(i,j) .gt. 0.0)imod0 = 1
       enddo
      enddo
      write(*,*) 'imod0=', imod0
      cldtc = cldtcm
c   Normalize re to ISCCP avg = (10.0+11.8)/2um = mean of Kiehl and Han
      sumre = 0.0
      nre = 0
      do j = 1,jcd
       do i = 1,icd
        if (cdrem(i,j) .ge. 0.0)then
         sumre = sumre + cdrem(i,j)
         nre = nre+1
        endif
       enddo
      enddo
      avgre = sumre/float(nre)
      avgisccp = (10.0+11.8)/2.0
      refac = avgisccp/avgre
      do j = 1,jm
       do i = 1,im
        if (cdrem(i,j) .ge. 0.0)then
         cdre(i,j) = refac*cdrem(i,j)
        else
         cdre(i,j) = cdrem(i,j)
        endif
       enddo
      enddo
c  File name: ISCCP clouds -- use preferentially for cover and LWP 
      call readnc_global_float_2d('clouds.nc','ccov', ccovi)
      call readnc_global_float_2d('clouds.nc','rlwp', rlwpi)

      ccov = ccovi
      rlwp = rlwpi
      if (ierrmodcld .eq. 1)then
       cdre = -1.0E10
      endif
c
      covmin = 9999.0
      covmax = -9999.0
      taumin = 9999.0
      taumax = -9999.0
      rlwmin = 9999.0
      rlwmax = -9999.0
      cremin = 9999.0
      cremax = -9999.0
      rdymin = 9999.0
      rdymax = -9999.0
      do jc = 1,jcd
       do ic = 1,icd
        rdays = float(ndmo)*ccov(ic,jc)*0.01
        ndy(ic,jc) = nint(rdays)
        if (ccov(ic,jc) .lt. 0.0)then
         jc1 = max(jc-1,1)
         ccov(ic,jc) = ccov(ic,jc1)
        endif
        if (ccov(ic,jc) .lt. 0.0)then
         jc1 = jc+1
         jc1 = min(jc1,jcd)
         do while (ccov(ic,jc1) .lt. 0.0 .and. jc1 .lt. jcd)
          jc1 = jc1+1
         enddo
         ccov(ic,jc) = ccov(ic,jc1)
        endif
        if (cldtc(ic,jc) .lt. 0.0)then
         jc1 = max(jc-1,1)
         cldtc(ic,jc) = cldtc(ic,jc1)
        endif
        if (cldtc(ic,jc) .lt. 0.0)then
         jc1 = jc+1
         jc1 = min(jc1,jcd)
         do while (cldtc(ic,jc1) .lt. 0.0 .and. jc1 .lt. jcd)
          jc1 = jc1+1
         enddo
         cldtc(ic,jc) = cldtc(ic,jc1)
        endif
        if (rlwp(ic,jc) .lt. 0.0)then
         jc1 = max(jc-1,1)
         rlwp(ic,jc) = rlwp(ic,jc1)
        endif
        if (rlwp(ic,jc) .lt. 0.0)then
         jc1 = jc+1
         jc1 = min(jc1,jcd)
         do while (rlwp(ic,jc1) .lt. 0.0 .and. jc1 .lt. jcd)
          jc1 = jc1+1
         enddo
         rlwp(ic,jc) = rlwp(ic,jc1)
        endif
        if (ierrmodcld .eq. 0)then
         if (cdre(ic,jc) .lt. 0.0)then
          jc1 = max(jc-1,1)
          cdre(ic,jc) = cdre(ic,jc1)
         endif
         if (cdre(ic,jc) .lt. 0.0)then
          jc1 = jc+1
          jc1 = min(jc1,jcd)
          do while (cdre(ic,jc1) .lt. 0.0 .and. jc1 .lt. jcd)
           jc1 = jc1+1
          enddo
          cdre(ic,jc) = cdre(ic,jc1)
         endif
        endif
c   find min/max
        if (ccov(ic,jc) .lt. covmin)then
         icmin = ic
         jcmin = jc
         covmin = ccov(ic,jc)
        endif
        if (ccov(ic,jc) .gt. covmax)then
         icmax = ic
         jcmax = jc
         covmax = ccov(ic,jc)
        endif
        if (cldtc(ic,jc) .lt. taumin)then
         itmin = ic
         jtmin = jc
         taumin = cldtc(ic,jc)
        endif
        if (cldtc(ic,jc) .gt. taumax)then
         itmax = ic
         jtmax = jc
         taumax = cldtc(ic,jc)
        endif
        if (rlwp(ic,jc) .lt. rlwmin)then
         ipmin = ic
         jpmin = jc
         rlwmin = rlwp(ic,jc)
        endif
        if (cdre(ic,jc) .gt. cremax)then
         iemax = ic
         jemax = jc
         cremax = cdre(ic,jc)
        endif
        if (cdre(ic,jc) .lt. cremin)then
         iemin = ic
         jemin = jc
         cremin = cdre(ic,jc)
        endif
        if (rlwp(ic,jc) .gt. rlwmax)then
         ipmax = ic
         jpmax = jc
         rlwmax = rlwp(ic,jc)
        endif
        if (ndy(ic,jc) .lt. rdymin)then
         inmin = ic
         jnmin = jc
         rdymin = ndy(ic,jc)
        endif
        if (ndy(ic,jc) .gt. rdymax)then
         inmax = ic
         jnmax = jc
         rdymax = ndy(ic,jc)
        endif
       enddo
      enddo
      if (taumax .lt. 0.0)then
       cdre = -1.0E10
       ierrmodcld = 1
      endif
c
      write(6,*)'Min/Max for Month ',imo
      write(6,10)'ccov  = ',covmin,icmin,jcmin,covmax,icmax,jcmax
      write(6,10)'cldtc = ',taumin,itmin,jtmin,taumax,itmax,jtmax
      write(6,10)'rlwp  = ',rlwmin,ipmin,jpmin,rlwmax,ipmax,jpmax
      if (ierrmodcld .eq. 0)then
       write(6,10)'cdre  = ',cremin,iemin,jemin,cremax,iemax,jemax
      endif
      write(6,10)'ndy   = ',rdymin,inmin,jnmin,rdymax,inmax,jnmax
      if (covmax .le. 0.0)then
       write(6,*)'Invalid cloud data'
       stop
      endif
c
c  Read in cloudy day distributions
      open(4,file='bcs/clddays.dat',status='old',form='formatted')
      do imt = 28,ndmo
       do icdy = 1,imt
        read(4,'(31i3)')(icld(icdy,it),it=1,31)
c        write(6,'(i3,31i3)')icdy,(icld(icdy,it),it=1,31)
       enddo
      enddo
      close(4)
c
c
c *****************************************************************
c  Read in other atmospheric optical parameters
      call readnc_global_float_2d('opt.nc','slporg', slporg)
      call readnc_global_float_2d('opt.nc','wsmorg', wsmorg)
      call readnc_global_float_2d('opt.nc','rhorg', rhorg)
      call readnc_global_float_2d('opt.nc','ozorg', ozorg)
      call readnc_global_float_2d('opt.nc','wvorg', wvorg)
c     open(4,file='opt.dat',status='old',form='unformatted')
c     read(4)imc,jmc
c     do m = 1,imo
c      read(4) slporg,wsmorg,rhorg,ozorg,wvorg
c     enddo
c     close(4)
      am = 1.0
      Vi = 25.0
c
      slmax = -99999.0
      slmin =  99999.0
      wsmax = -99999.0
      wsmin =  99999.0
      rhmax = -99999.0
      rhmin =  99999.0
      ozmax = -99999.0
      ozmin =  99999.0
      wvmax = -99999.0
      wvmin =  99999.0
      do j = 1,jatm
       do i = 1,iatm
        slp(i,j) = slporg(i,j)
        if (slp(i,j) .gt. flag)then
         if (slp(i,j) .gt. slmax)then
          islmax = i
          jslmax = j
          slmax = slp(i,j)
         endif
         if (slp(i,j) .lt. slmin)then
          islmin = i
          jslmin = j
          slmin = slp(i,j)
         endif
        endif
        if (wsmorg(i,j) .gt. flag)then
         if (wsmorg(i,j) .gt. wsmax)then
          iwsmax = i
          jwsmax = j
          wsmax = wsmorg(i,j)
         endif
         if (wsmorg(i,j) .lt. wsmin)then
          iwsmin = i
          jwsmin = j
          wsmin = wsmorg(i,j)
         endif
        endif
        rhorg(i,j) = min(rhorg(i,j),100.0)  !prevent excess saturation
        if (rhorg(i,j) .gt. flag)then
         if (rhorg(i,j) .gt. rhmax)then
          irhmax = i
          jrhmax = j
          rhmax = rhorg(i,j)
         endif
         if (rhorg(i,j) .lt. rhmin)then
          irhmin = i
          jrhmin = j
          rhmin = rhorg(i,j)
         endif
        endif
        if (ozorg(i,j) .lt. 0.0)then
         j1 = max(j-1,1)
         ozorg(i,j) = ozorg(i,j1)
        endif
        if (ozorg(i,j) .lt. 0.0)then
         j1 = j+1
         j1 = min(j1,jatm)
         do while (ozorg(i,j1) .lt. 0.0 .and. j1 .lt. jatm)
          j1 = j1+1
         enddo
         ozorg(i,j) = ozorg(i,j1)
        endif
        if (ozorg(i,j) .gt. ozmax)then
         iozmax = i
         jozmax = j
         ozmax = ozorg(i,j)
        endif
        if (ozorg(i,j) .lt. ozmin)then
         iozmin = i
         jozmin = j
         ozmin = ozorg(i,j)
        endif
        if (wvorg(i,j) .gt. flag)then
         if (wvorg(i,j) .gt. wvmax)then
          iwvmax = i
          jwvmax = j
          wvmax = wvorg(i,j)
         endif
         if (wvorg(i,j) .lt. wvmin)then
          iwvmin = i
          jwvmin = j
          wvmin = wvorg(i,j)
         endif
        endif
        wsm(i,j) = wsmorg(i,j)
        rh(i,j) = rhorg(i,j)
        oz(i,j) = ozorg(i,j)
        wv(i,j) = wvorg(i,j)
       enddo
      enddo
c
      write(6,10)'slp  = ',slmin,islmin,jslmax,slmax,islmax,jslmax
      write(6,10)'wsm  = ',wsmin,iwsmin,jwsmin,wsmax,iwsmax,jwsmax
      write(6,10)'rh   = ',rhmin,irhmin,jrhmin,rhmax,irhmax,jrhmax
      write(6,10)'oz   = ',ozmin,iozmin,jozmin,ozmax,iozmax,jozmax
      write(6,10)'wv   = ',wvmin,iwvmin,jwvmin,wvmax,iwvmax,jwvmax
c
c******************************************************************
c  Read in aerosols
      taua = flagam
c   Use MODIS aerosols if available
      open(4,file='modaer.dat',status='old',form='unformatted',
     *     err=7)
      read(4)taua,asymp,ssalb
      close(4)
      go to 9
7     write(6,*)'modaer not found; use AVHRR aerosol tau'
c   Else use AVHRR aerosol tau and MODIS climatologies for spectal
c   variability of taua, asymmetry parameter (asymp) and single-
c   scattering albedo (ssalb)
      open(4,file='aot.dat',status='old',form='unformatted',
     *     err=8)
      do immo = 1,imo
       read(4)avhrrtau     !taua at 0.63um
      enddo
      close(4)
c   get monthly climatology for missing data
      open(4,file='aotclim.dat',status='old',form='unformatted',
     *     err=8)
      do immo = 1,imo
       read(4)avhrrtauclim     
      enddo
      close(4)
      do j = 1,jatm
       do i = 1,iatm
        if (avhrrtau(i,j) .lt. 0.0)then
         avhrrtau(i,j) = avhrrtauclim(i,j)
        endif
       enddo
      enddo
      write(6,*)'avhrrtau = ',avhrrtau(50,100)
      open(4,file='modaerclim.dat',status='old',form='unformatted',
     *     err=8)
      read(4)taua,asymp,ssalb
      close(4)
c   Closest OASIM band is 625 nm; AVHRR is 630 nm
      do j = 1,jatm
       do i = 1,iatm
        if (avhrrtau(i,j) .ge. 0.0)then
         do nl = 1,nlt
          facaer(nl) = taua(i,j,nl)/taua(i,j,14)
         enddo
c    subtract AVHRR bias (add a negative to a negative bias)
         taua(i,j,14) = avhrrtau(i,j) + 0.025   !625 nm
         do nl = 1,nlt
          taua(i,j,nl) = facaer(nl)*taua(i,j,14)
         enddo
        endif
       enddo
      enddo
8     continue
      close(4)
9     continue
c
      tamin =  99999.0
      tamax = -99999.0
      asmin =  99999.0
      asmax = -99999.0
      wamin =  99999.0
      wamax = -99999.0
      do j = 1,jatm
       do i = 1,iatm
        do nl = 1,nlt
         if (taua(i,j,nl) .lt. tamin)then
          itamin = i
          jtamin = j
          tamin = taua(i,j,nl)
         endif
         if (taua(i,j,nl) .gt. tamax)then
          itamax = i
          jtamax = j
          tamax = taua(i,j,nl)
         endif
         if (asymp(i,j,nl) .lt. asmin)then
          iasmin = i
          jasmin = j
          asmin = asymp(i,j,nl)
         endif
         if (asymp(i,j,nl) .gt. asmax)then
          iasmax = i
          jasmax = j
          asmax = asymp(i,j,nl)
         endif
         if (ssalb(i,j,nl) .lt. wamin)then
          iwamin = i
          jwamin = j
          wamin = ssalb(i,j,nl)
         endif
         if (ssalb(i,j,nl) .gt. wamax)then
          iwamax = i
          jwamax = j
          wamax = ssalb(i,j,nl)
         endif
        enddo
       enddo
      enddo
c
      write(6,10)'taua  = ',tamin,itamin,jtamax,tamax,itamax,jtamax
      write(6,10)'asymp = ',asmin,iasmin,jasmin,asmax,iasmax,jasmax
      write(6,10)'ssalb = ',wamin,iwamin,jwamin,wamax,iwamax,jwamax
c
10    format(a10,f8.2,1x,'(',i3,',',i3,')',3x,f8.2,1x,'(',i3,',',i3,')')
      return
      end
