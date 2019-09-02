      subroutine monavgpar(rod,ros,Ed,Es)
c
c  Computes monthly mean PAR. BELOW surface (0-).
c
#include "define.h"
#include "comlte.h"
      byte msk
      real rod(im,jm,nlt,nhn),ros(im,jm,nlt,nhn)
      real Ed(im,jm,nlt,nhn),Es(im,jm,nlt,nhn)
      real avgparq(im,jm)
      real area(jm)
      common /bwq/ WtoQ(nlt)
      common /bmsk/ msk(im,jm)
      data flagl /999.0E9/
c
      avgparq = 0.0
      rminq = 99999.0
      rmaxq = -99999.0
      Edq = flagl
      Esq = flagl
      do j = 1,jm
       do i = 1,im
        if (msk(i,j) .lt. 2)then
         do ihs = 1,nhn
          do nl = npst,npnd   !PAR range index 3-17
           Ed0m = Ed(i,j,nl,ihs)*(1.0-rod(i,j,nl,ihs))
           Es0m = Es(i,j,nl,ihs)*(1.0-ros(i,j,nl,ihs))
           Edq = Ed0m*WtoQ(nl)*2.0*3600.0    !W/m2 to quanta/m2/hr
           Esq = Es0m*WtoQ(nl)*2.0*3600.0    !W/m2 to quanta/m2/hr
           totq = Edq + Esq
           avgparq(i,j) = avgparq(i,j) + totq
          enddo
         enddo
        else
         avgparq(i,j) = flagl
        endif
       enddo
      enddo
      call getarea(area)
      gloavgparq = 0.0
      gloarea = 0.0
      do j = 1,jm
       do i = 1,im
        if (msk(i,j) .lt. 2)then
         gloavgparq = gloavgparq + avgparq(i,j)*area(j)
         gloarea = gloarea + area(j)
         rminq = min(rminq,avgparq(i,j))
         rmaxq = max(rmaxq,avgparq(i,j))
        endif
       enddo
      enddo
      gloavgparq = gloavgparq/gloarea
      write(6,*)'Mean,min,max PAR moles quanta/m2/day = ',gloavgparq,
     *rminq,rmaxq
c
      open(3,file='avgparq.dat',status='unknown',form='unformatted')
      write(3)avgparq
      close(3)
c
      call WRITE_AVE_2D_1v('parq.nc', 'PARQ', avgparq)
      end
