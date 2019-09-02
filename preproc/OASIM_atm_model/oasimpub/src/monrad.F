      program monrad
c
c  Computes mean monthly spectral irradiance.
c
#include "define.h"
#include "rad.h"
#include "comlte.h"
#include "comdata.h"
      real avgd(im,jm,nlt,nhn),avgs(im,jm,nlt,nhn)
      real avgd_0m(im,jm,nlt,nhn),avgs_0m(im,jm,nlt,nhn)
      integer*2 iavgd(im,jm,nlt,nhn),iavgs(im,jm,nlt,nhn)
      real solz(im,jm)
      real rod(im,jm,nlt),ros(im,jm,nlt)
      real rodihs(im,jm,nlt,nhn),rosihs(im,jm,nlt,nhn)
      integer nihs(im,jm,nlt,nhn)
      data rminval /1.0E-36/
c
      open(4,file='day.dat',status='old',form='formatted')
      read(4,*) iday_bioptimod
      write(*,*) iday_bioptimod
      close(4)
      open(4,file='month.dat',status='old',form='formatted')
      read(4,*) imo
      close(4)
      open(4,file='year.dat',status='old',form='formatted')
      read(4,*) iyr
      close(4)
      idst = iday_bioptimod
      iend = iday_bioptimod
      
      call setdate(idst,imo,iyr,istday,ndmo,lpyr)
      write(6,*)'Start day, month, year = ',idst,imo,iyr

#if RAD2HR
      write(*,*) 'RAD2HR Active'
#endif
c  BIOLOGICAL/LIGHT SETUP
c  Read in constants, light data
        call setlte(1,ndmo)
c
      ndmo  = 1
      write(6,*)'ndmo = ',ndmo, 'only one day'
 
c  Initialize 
      avgd = 0.0
      avgs = 0.0
      avgd_0m = 0.0
      avgs_0m = 0.0
      iavgd = 0
      iavgs = 0
      rodihs = 0.0
      rosihs = 0.0
      nihs = 0
c
      npt = 5
      jt = jm/2
      it = im/2
      iht = 7
!-----------------  MAIN LOOP -------------------------------------
      do idmo = iday_bioptimod,iday_bioptimod
       iday = idmo
cbioptimod  iday = istday + (idmo-1)
       write(6,*)'idmo,iday = ',idmo,iday

c  Hour loop
       do ihr = 1,24
        write(6,*)'idmo,ihr = ',idmo,ihr
        ihtst = ihr+1
#if RAD2HR
        if (mod(ihtst,2) .eq. 0)then
         ihs = int(ihtst/2)
         call sfcsolz(iyr,iday,ihr,solz)
         do j = 1,jm
          do i = 1,im
           call ocalbedo(rad,solz(i,j),wsm(i,j),rod(i,j,:),ros(i,j,:))
           do nl = 1,nlt
            rodihs(i,j,nl,ihs) = rodihs(i,j,nl,ihs) + rod(i,j,nl)
            rosihs(i,j,nl,ihs) = rosihs(i,j,nl,ihs) + ros(i,j,nl)
            nihs(i,j,nl,ihs) = nihs(i,j,nl,ihs)+1
           enddo
          enddo
         enddo
         call sfcirr(idmo,iyr,iday,ihr,imo,solz)
#else
         ihs = ihr
         call sfcsolz(iyr,iday,ihr,solz)
         do j = 1,jm
          do i = 1,im
           call ocalbedo(rad,solz(i,j),wsm(i,j),rod(i,j,:),ros(i,j,:))
           do nl = 1,nlt
            rodihs(i,j,nl,ihs) = rodihs(i,j,nl,ihs) + rod(i,j,nl)
            rosihs(i,j,nl,ihs) = rosihs(i,j,nl,ihs) + ros(i,j,nl)
            nihs(i,j,nl,ihs) = nihs(i,j,nl,ihs)+1
           enddo
          enddo
         enddo
         call sfcirr(idmo,iyr,iday,ihr,imo,solz)
#endif
c   Sum as go along
         do nl = 1,nlt
          do j = 1,jm
           do i = 1,im
            avgd(i,j,nl,ihs) = avgd(i,j,nl,ihs) + Eda(i,j,nl)
            avgs(i,j,nl,ihs) = avgs(i,j,nl,ihs) + Esa(i,j,nl)
            avgd_0m(i,j,nl,ihs) = avgd_0m(i,j,nl,ihs) + Eda(i,j,nl)
     &                            *( 1.0 -rod(i,j,nl) ) 
            avgs_0m(i,j,nl,ihs) = avgs_0m(i,j,nl,ihs) + Esa(i,j,nl)
     &                            *( 1.0 -ros(i,j,nl) ) 
           enddo
          enddo
         enddo
       write(6,*)'ihr,it,npt,Eda,Esa = ',ihr,it,npt,Eda(it,jt,npt),
     *Esa(it,jt,npt)
       write(6,*)'ihr,it,npt,avgd,avgs = ',ihr,it,npt,
     *avgd(it,jt,npt,ihs),avgs(it,jt,npt,ihs)
#if RAD2HR
        endif
#endif
       enddo    !end of hour loop
c
      enddo   !end of day loop
c
c  Take mean over each hour for the month
      write(*,*) 'ndmo', ndmo
      do ihs = 1,nhn
       do nl = 1,nlt
        do j = 1,jm
         do i = 1,im
          avgd(i,j,nl,ihs) = avgd(i,j,nl,ihs)/float(ndmo)
          avgd(i,j,nl,ihs) = max(avgd(i,j,nl,ihs),rminval)
          avgs(i,j,nl,ihs) = avgs(i,j,nl,ihs)/float(ndmo)
          avgs(i,j,nl,ihs) = max(avgs(i,j,nl,ihs),rminval)

          avgd_0m(i,j,nl,ihs) = avgd_0m(i,j,nl,ihs)/float(ndmo)
          avgd_0m(i,j,nl,ihs) = max(avgd_0m(i,j,nl,ihs),rminval)
          avgs_0m(i,j,nl,ihs) = avgs_0m(i,j,nl,ihs)/float(ndmo)
          avgs_0m(i,j,nl,ihs) = max(avgs_0m(i,j,nl,ihs),rminval)

          iavgd(i,j,nl,ihs) = nint(avgd(i,j,nl,ihs)*200.0)
          iavgs(i,j,nl,ihs) = nint(avgs(i,j,nl,ihs)*200.0)
          if (iavgd(i,j,nl,ihs) .gt. 32700)then
           write(6,*)'ERROR -- scale factor too big'
           write(6,*)'i,j,nl,ihs,avgd,iavgd = ',i,j,nl,ihs,
     *      avgd(i,j,nl,ihs),iavgd(i,j,nl,ihs)
          endif
          if (iavgs(i,j,nl,ihs) .gt. 32700)then
           write(6,*)'ERROR -- scale factor too big'
           write(6,*)'i,j,nl,ihs,avgs,iavgs = ',i,j,nl,ihs,
     *      avgs(i,j,nl,ihs),iavgs(i,j,nl,ihs)
          endif
         enddo
        enddo
       enddo
      enddo
      write(6,*)'hour,lam,avgd,avgs,avgd_0m,avgs_0m'
      do ih = 1,nhn
      write(6,*)ih,lam(npt),avgd(it,jt,npt,ih),
     *avgs(it,jt,npt,ih),avgd_0m(it,jt,npt,ih),
     *avgs_0m(it,jt,npt,ih)
      enddo
      totqt = 0.0
      do nl = 1,nlt
       totqt = totqt + avgd(it,jt,nl,iht)
     *               + avgs(it,jt,nl,iht)
      enddo
      write(6,*)'Average total irrad for this hour = ',totqt
c
c  Write to data file
      open(4,file='rad',status='unknown',form='unformatted')
      write(4)iavgd,iavgs
      close(4)
      call  WRITE_AVE_2D_2v_l_h('rad_0p.nc', 'Ed_0p', 'Es_0p'
     &                         , avgd, avgs)    
      call  WRITE_AVE_2D_2v_l_h('rad_0m.nc', 'Ed_0m', 'Es_0m'
     &                         , avgd_0m, avgs_0m)    
c
c   Mean surface reflectance
      do ihs = 1,nhn
       do j = 1,jm
        do i = 1,im
         do nl = 1,nlt
          denom = float(nihs(i,j,nl,ihs))
          rodihs(i,j,nl,ihs) = rodihs(i,j,nl,ihs)/denom
          rosihs(i,j,nl,ihs) = rosihs(i,j,nl,ihs)/denom
         enddo
        enddo
       enddo
      enddo
      call monavgswr(rodihs,rosihs,avgd,avgs)
      call monavgpar(rodihs,rosihs,avgd,avgs)
      call monavgday0m(rodihs,rosihs,avgd,avgs)
c
      end
