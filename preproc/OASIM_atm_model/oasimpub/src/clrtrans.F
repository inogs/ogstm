      subroutine clrtrans(cosunz,rm,rmp,rmo,ws,relhum,am,Vi,
     *ta,wa,asym,Td,Ts)
c
c  Model for atmospheric transmittance of solar irradiance through
c  a cloudless maritime atmosphere.  Computes direct and diffuse 
c  separately.  From Gregg and Carder (1990) Limnology and 
c  Oceanography 35(8): 1657-1675.
c
c  Td is spectral clear sky direct transmittance
c  Ts is spectral clear sky diffuse transmittance
c
#include "comlte.h"
      save ifst,rlamu
      real ta(nlt),asym(nlt),wa(nlt)
      real rlamu(nlt)
      real Td(nlt),Ts(nlt)
      data ifst /0/
c
c  First time only
      if (ifst .eq. 0)then
       do nl = 1,nlt
        rlamu(nl) = float(lam(nl))*1.0E-3    !lambda in um
        Td(nl) = 0.0
        Ts(nl) = 0.0
       enddo
       ifst = 1
      endif
c
c  Obtain aerosol parameters; simplified Navy aerosol model
      call navaer(relhum,am,Vi,ws,beta,eta,wa1,afs,bfs)
c
c  Compute spectral transmittance
      do nl = 1,nlt
c    Rayleigh
       rtra = exp(-thray(nl)*rmp)       !transmittance
c   Aerosols
       if (ta(nl) .lt. 0.0)then
        ta(nl) = beta*rlamu(nl)**eta
       endif
       if (wa(nl) .lt. 0.0)then
        omegaa = wa1
       else
        omegaa = wa(nl)
       endif
       if (asym(nl) .ge. 0.0)then
c bioptimod added  + 0.01 to avoid NaN
        alg = log(1.0-asym(nl)+0.01)
        afs = alg*(1.459+alg*(.1595+alg*.4129))
        bfs = alg*(.0783+alg*(-.3824-alg*.5874))
       endif
       if (ta(nl) .lt. 0.0 .or. omegaa .lt. 0.0)then
        write(6,*)'ERROR in ta or omegaa'
        write(6,*)'nl,ta,wa,asym = ', nl,ta(nl),wa(nl),asym(nl)
       endif
       Fa = 1.0 - 0.5*exp((afs+bfs*cosunz)*cosunz)
       if (Fa .lt. 0.0)then
        write(6,*)'ERROR in Fa'
        write(6,*)'nl,ta,wa,asym = ', nl,ta(nl),wa(nl),asym(nl)
       endif
       tarm = ta(nl)*rm
       atra = exp(-tarm)
       taa = exp(-(1.0-omegaa)*tarm)
       tas = exp(-omegaa*tarm)
c  Direct transmittance
       Td(nl) = rtra*atra
c
c   Diffuse transmittance
       dray = taa*0.5*(1.0-rtra**.95)
       daer = rtra**1.5*taa*Fa*(1.0-tas)
c
c  Total diffuse
       Ts(nl) = dray + daer
c
      enddo
c
      return
      end
c
c *****************************************************************
      subroutine navaer(relhum,am,Vi,ws,beta,eta,wa,afs,bfs)
c
c  Computes aerosol parameters according to a simplified version
c  of the Navy marine aerosol model.
c
      real a(3),ro(3),dndr(3),r(3)
      data ro /0.03,0.24,2.0/
      data r /0.1,1.0,10.0/
      data rlam /0.55/
c
      wsm = ws
c
c  Relative humidity factor
c      if (relhum .ge. 100.0)relhum = 99.9
      relhum = min(99.9,relhum)
      rnum = 2.0 - relhum/100.0
      rden = 6.0*(1.0-relhum/100.0)
      frh = (rnum/rden)**0.333
c
c  Size distribution amplitude components
      a(1) = 2000.0*am*am
      a(2) = 5.866*(wsm-2.2)
c      if (a(2) .lt. 0.5)a(2) = 0.5
      a(2) = max(0.5,a(2))
      a(3) = 0.01527*(ws-2.2)*0.05        !from Hughes 1987
c      if (a(3) .lt. 1.4E-5)a(3) = 1.4E-5
      a(3) = max(1.4E-5,a(3))
c
c  Compute size distribution at three selected radii according to
c  Navy method
      do n = 1,3
       dndr(n) = 0.0
       do i = 1,3
        rden = frh*ro(i)
        arg = log(r(n)/rden)*log(r(n)/rden)
        rval = a(i)*exp(-arg)/frh
        dndr(n) = dndr(n) + rval
       enddo
      enddo
c
c  Least squares approximation
      sumx = 0.0
      sumy = 0.0
      sumxy = 0.0
      sumx2 = 0.0
      do n = 1,3
       rlrn = log10(r(n))
       rldndr = log10(dndr(n))
       sumx = sumx + rlrn
       sumy = sumy + rldndr
       sumxy = sumxy + rlrn*rldndr
       sumx2 = sumx2 + rlrn*rlrn
      enddo
      gama = sumxy/sumx2
      rlogc = sumy/3.0 - gama*sumx/3.0
      alpha = -(gama+3.0)
      eta = -alpha
c
c  Compute beta
      cext = 3.91/Vi
      beta = cext*rlam**alpha
c
c  Compute asymmetry parameter -- a function of alpha
      if (alpha .gt. 1.2)then
       asymp = 0.65
      else if (alpha .lt. 0.0)then
       asymp = 0.82
      else
       asymp = -0.14167*alpha + 0.82
      endif
c
c  Forward scattering coefficients
      alg = log(1.0-asymp)
      afs = alg*(1.459+alg*(.1595+alg*.4129))
      bfs = alg*(.0783+alg*(-.3824-alg*.5874))
c
c  Single scattering albedo at 550; function of RH
      wa = (-0.0032*am + 0.972)*exp(3.06E-4*relhum)
c
      return
      end
c
c *****************************************************************
      subroutine sfcrfl(rad,theta,ws,rod,ros)
c
c  Computes surface reflectance for direct (rod) and diffuse (ros)
c  components separately, as a function of theta, wind speed or
c  stress.
c  Includes spectral dependence of foam reflectance derived from Frouin
c  et al., 1996 (JGR)
c
#include "comlte.h"
      save ifst,rn,roair
      real rod(nlt),ros(nlt)
      common /bwfac/ wfac(nlt)
      data a0,a1,a2,a3 /0.9976,0.2194,5.554E-2,6.7E-3/
      data b0,b1,b2,b3 /5.026,-0.01138,9.552E-6,-2.698E-9/
      data ifst /0/
c
      if (ifst .eq. 0)then
       rn = 1.341    !index of refraction of pure seawater
       roair = 1.2E3     !density of air g/m3
       do nl = 1,nlt
        rlam = float(lam(nl))
        if (lam(nl) .lt. 900)then
         t = exp(-(aw(nl)+0.5*bw(nl)))
         tlog = alog(1.0E-36+t)
         fac = a0 + a1*tlog + a2*tlog*tlog + a3*tlog*tlog*tlog
         wfac(nl) = min(fac,1.0)
         wfac(nl) = max(wfac(nl),0.0) !BIOPTIMOD corection: from fac to wfac
        else
         fac = b0 + b1*rlam + b2*rlam*rlam + b3*rlam*rlam*rlam
         wfac(nl) = max(fac,0.0)
        endif
       enddo
       ifst = 1
      endif
c
c  Foam and diffuse reflectance
      if (ws .gt. 4.0)then
       if (ws .le. 7.0)then
        cn = 6.2E-4 + 1.56E-3/ws
        rof = roair*cn*2.2E-5*ws*ws - 4.0E-4
       else
        cn = 0.49E-3 + 0.065E-3*ws
        rof = (roair*cn*4.5E-5 - 4.0E-5)*ws*ws
       endif
       rosps = 0.057
      else
       rof = 0.0
       rosps = 0.066
      endif
c
c  Direct
c   Fresnel reflectance for theta < 40, ws < 2 m/s
      if (theta .lt. 40.0 .or. ws .lt. 2.0)then
       if (theta .eq. 0.0)then
        rospd = 0.0211
       else
        rtheta = theta/rad
        sintr = sin(rtheta)/rn
        rthetar = asin(sintr)
        rmin = rtheta - rthetar
        rpls = rtheta + rthetar
        sinrmin = sin(rmin)
        sinrpls = sin(rpls)
        tanrmin = tan(rmin)
        tanrpls = tan(rpls)
        sinp = (sinrmin*sinrmin)/(sinrpls*sinrpls)
        tanp = (tanrmin*tanrmin)/(tanrpls*tanrpls)
        rospd = 0.5*(sinp + tanp)
       endif
      else
c   Empirical fit otherwise
       a = 0.0253
       b = -7.14E-4*ws + 0.0618
       rospd = a*exp(b*(theta-40.0))
      endif
c
c  Reflectance totals
      do nl = 1,nlt
       rod(nl) = rospd + rof*wfac(nl)
       ros(nl) = rosps + rof*wfac(nl)
      enddo
c
      return
      end
