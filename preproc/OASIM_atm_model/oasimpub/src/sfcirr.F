      subroutine sfcirr(idmo,iyr,iday,ihr,imo,solz)
c
c  Computes hourly integrated iradiance given location and time.
c
#include "define.h"
#include "rad.h"
#include "comlte.h"
#include "comdata.h"
      save iprt,o24,ifst
      real solz(im,jm)
      real ta(nlt),wa(nlt),asym(nlt)
      real Ed(nlt),Es(nlt)
      data ifst /0/
c
c  Set up
      if (ifst .eq. 0)then
#if RAD2HR
c  compute irradiance every 2 hours
       iprt = 13
#else
       iprt = 12
#endif
       o24 = 1.0/24.0
       ifst = 1
      endif
c
c  Initialize arrays everytime
      Eda = 0.0
      Esa = 0.0
c
c  Obtain Earth-Sun distance only for the middle of the hour
      hr = float(ihr)
      rday = float(iday) + hr*o24
      daycor = (1.0+1.67E-2*cos(pi2*(rday-3.0)/365.0))**2
c
c  Compute Ed, Es just above surface
      do jl = 1,jm
        do il = 1,im
           cosunz = cos(solz(il,jl)/rad)
           sunz = solz(il,jl)
           if (sunz .lt. 90.0)then
c   atmospheric gases
              pres   = slp(il,jl)
              ws     = wsm(il,jl)
              ozone  = oz(il,jl)
              wvapor = wv(il,jl)
              relhum = rh(il,jl)
c   aerosols
              do nl = 1,nlt
                 ta(nl) = taua(il,jl,nl)
                 asym(nl) = asymp(il,jl,nl)
                 wa(nl) = ssalb(il,jl,nl)
              enddo
c   clouds
              call findxy(im,jm,il,jl,xl,yl)
              call findbin(icd,jcd,xl,yl,ic,jc)
              cov = ccov(ic,jc)
              cldtau = cldtc(ic,jc)
              clwp = rlwp(ic,jc)
              re = cdre(ic,jc)
c   Ed, Es
              call light(sunz,cosunz,daycor,pres,ws,ozone,wvapor,relhum,
     *        ta,wa,asym,am,Vi,cov,cldtau,clwp,re,
     *        Ed,Es)
              sirr = 0.0
              do nl = 1,nlt
                 Eda(il,jl,nl) = Ed(nl)
                 Esa(il,jl,nl) = Es(nl)
                 sirr = sirr + Eda(il,jl,nl) + Esa(il,jl,nl)
              enddo
        else
          sirr = 0.0
          Eda(il,jl,:) = 0.0
          Esa(il,jl,:) = 0.0
        endif
c       write(6,*)'Surface irrad W/m2 = ',jl,sirr
        enddo   !end of il loop
      enddo    !end of jl loop
c
c  Noon print-out
c      if (ihr .eq. iprt)then
       jl = jm/2
       il = im/2
        sirr = 0.0
        do nl = 1,nlt
         sirr = sirr + Eda(il,jl,nl) + Esa(il,jl,nl)
        enddo
        write(6,'(a30,4i6,f10.2)')'Surface irrad W/m2 = ',
     *idmo,ihr,jl,il,sirr
c       enddo
c      enddo
c      endif
c
      return
      end
