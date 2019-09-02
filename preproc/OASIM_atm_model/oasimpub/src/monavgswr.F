      subroutine monavgswr(rod,ros,Ed,Es)
c
c  Computes monthly mean irradiances or total shortwave radiation.
c  Accounts for surface reflectance.
c
#include "define.h"
#include "comlte.h"
      byte msk
      real Ed(im,jm,nlt,nhn),Es(im,jm,nlt,nhn)
      real rod(im,jm,nlt,nhn),ros(im,jm,nlt,nhn)
      real avgswr(im,jm)
      real area(jm)
      common /bmsk/ msk(im,jm)
      data flagl /999.0E9/
c
c  mask.dat: 1-degree, 0=open ocean, 1=coast, 2=land
      open(4,file='mask.dat',status='old',form='unformatted')
      read(4)msk
      close(4)
c
      avgswr = 0.0
      do j = 1,jm
       do i = 1,im
        do ihs = 1,nhn
         do nl = 1,nlt
          Ed0m = 0.0
          Es0m = 0.0
          Ed0m = Ed(i,j,nl,ihs)*(1.0-rod(i,j,nl,ihs))
          Es0m = Es(i,j,nl,ihs)*(1.0-ros(i,j,nl,ihs))
          tot = Ed0m + Es0m
          avgswr(i,j) = avgswr(i,j) + tot
         enddo
        enddo
       enddo
      enddo
      rmin = 99999.0
      rmax = -99999.0
      gloavg = 0.0
      gloarea = 0.0
      call getarea(area)
      do j = 1,jm
       do i = 1,im
        avgswr(i,j) = avgswr(i,j)/float(nhn)
        if (msk(i,j) .lt. 2)then
         gloavg = gloavg + avgswr(i,j)*area(j)
         gloarea = gloarea + area(j)
         rmin = min(rmin,avgswr(i,j))
         rmax = max(rmax,avgswr(i,j))
        endif
       enddo
      enddo
      gloavg = gloavg/gloarea
      write(6,*)'Mean,min,max Irrad W/m2 = ',gloavg,rmin,rmax
c
      open(3,file='avgirr.dat',status='unknown',form='unformatted')
      write(3)avgswr
      close(3)

      call WRITE_AVE_2D_1v('swr.nc', 'SWR', avgswr)

c
      return
      end
      
