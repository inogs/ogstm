      subroutine setlte(imo,ndmo)
c
c  Computes constants for global irradiance calculations, reads in 
c  required data files, and otherwise obtains one-time-only information 
c  necessary for the run.
c
#include "rad.h"
#include "comlte.h"
      common /bwq/ WtoQ(nlt)
c
c  Degrees to radians conversion
      pi = dacos(-1.0D0)
      pi2 = pi*2.0
      rad = 180.0D0/pi
c
c  Obtain Light data
      call lidata
      call rdatopt(imo,ndmo)
c
c  Quanta conversion
      h = 6.6256E-34   !Plancks constant J sec
      c = 2.998E8      !speed of light m/sec
      hc = 1.0/(h*c)
      oavo = 1.0/6.023E23   ! 1/Avogadros number
      hcoavo = hc*oavo
      do nl = 1,nlt
       rlamm = float(lam(nl))*1.0E-9  !lambda in m
       WtoQ(nl) = rlamm*hcoavo        !Watts to quanta conversion
      enddo
c
      return
      end
