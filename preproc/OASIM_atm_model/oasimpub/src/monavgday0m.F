      subroutine monavgday0m(rod,ros,Ed,Es)
c
c  Computes monthly mean non-diurnal spectral irradiance.
c  Accounts for surface reflectance.  
c
#include "define.h"
#include "comlte.h"
      real Ed(im,jm,nlt,nhn),Es(im,jm,nlt,nhn)
      real rod(im,jm,nlt,nhn),ros(im,jm,nlt,nhn)
      real Edday(im,jm,nlt),Esday(im,jm,nlt)
      real avgday(im,jm)
      real area(jm)
      byte msk
      common /bmsk/ msk(im,jm)
      data flagl /999.0E9/
c
      Edday = flagl
      Esday = flagl
      do j = 1,jm
       do i = 1,im
        if (msk(i,j) .lt. 2)then
         do nl = 1,nlt
          Ed0m = 0.0
          Es0m = 0.0
          do ihs = 1,nhn
           Ed0m = Ed0m + Ed(i,j,nl,ihs)*(1.0-rod(i,j,nl,ihs))
           Es0m = Es0m + Es(i,j,nl,ihs)*(1.0-ros(i,j,nl,ihs))
          enddo
          Edday(i,j,nl) = Ed0m/float(nhn)
          Esday(i,j,nl) = Es0m/float(nhn)
         enddo
        endif
       enddo
      enddo
      avgday = 0.0
      rmin = 99999.0
      rmax = -99999.0
      gloavg = 0.0
      gloarea = 0.0
      call getarea(area)
      do j = 1,jm
       do i = 1,im
        if (msk(i,j) .lt. 2)then
         do nl = 1,nlt
          avgday(i,j) = avgday(i,j) + Edday(i,j,nl) + Esday(i,j,nl)
         enddo
         gloavg = gloavg + avgday(i,j)*area(j)
         gloarea = gloarea + area(j)
         rmin = min(rmin,avgday(i,j))
         rmax = max(rmax,avgday(i,j))
        endif
       enddo
      enddo
      gloavg = gloavg/gloarea
      write(6,*)'Mean,min,max Irrad0m W/m2 = ',gloavg,rmin,rmax
c
      open(3,file='eds.dat',status='unknown',form='unformatted')
      write(3)Edday,Esday
      close(3)
c
      call WRITE_AVE_2D_2v_l('eds.nc', 'Ed_0m', 'Es_0m', Edday, Esday)
      return
      end
      
