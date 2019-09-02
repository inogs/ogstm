      subroutine light_verbose(sunz,cosunz,daycor,pres,ws,ozone,wvapor,relhum,
     * ta,wa,asym,am,Vi,cov,cldtau,clwp,re,
     * Ed,Es)
c
c  Calls clrtrans.F to get cloud-free transmittance and slingo.F to
c  get cloudy transmittance, then computes total irradiance in
c  W/m2/(variable)nm weighted by the cloudiness.
c
c  Tdclr is spectral clear sky direct transmittance
c  Tsclr is spectral clear sky diffuse transmittance
c  Tdcld is spectral cloudy direct transmittance
c  Tscld is spectral cloudy diffuse transmittance
c
#include "comlte.h"
      save ifst,p0,ozfac1,ozfac2
      parameter(ncld=24)
      real ta(nlt),wa(nlt),asym(nlt),Tgas(nlt)
      real Tdclr(nlt),Tsclr(nlt)
      real Tdcld(ncld),Tscld(ncld)
      real Ed(nlt),Es(nlt)
      real Edclr(nlt),Esclr(nlt)
      real Edcld(nlt),Escld(nlt)
      common /bica/ ica(nlt)
      data ifst /0/
c
c  Initial constants
      if (ifst .eq. 0)then
       ozfac1 = 44.0/6370.0
       ozfac2 = 1.0 + 22.0/6370.0
       p0 = 1013.25
       ifst = 1
      endif
c
      if (pres .lt. 0.0 .or. ws .lt. 0.0 .or. relhum .lt. 0.0
     *.or. ozone .lt. 0.0 .or. wvapor .lt. 0.0)then
       Ed = 0.0
       Es = 0.0
       go to 5
      endif
c
c  Compute atmospheric path lengths (air mass); not pressure-corrected
      rtmp = (93.885-sunz)**(-1.253)
      rmu0 = cosunz+0.15*rtmp
      rm = 1.0/rmu0
      otmp = (cosunz*cosunz+ozfac1)**0.5
      rmo = ozfac2/otmp
c
c  Compute pressure-corrected atmospheric path length (air mass)
      rmp = pres/p0*rm
c
c  Loop to compute total irradiances at each grid point
c   Compute direct and diffuse irradiance for a cloudy and non-cloudy
c   atmosphere
c   Account for gaseous absorption
      do nl = 1,nlt
c    Ozone
       to = oza(nl)*ozone*0.001   !optical thickness
       oarg = -to*rmo
c   Oxygen/gases
       ag = ao(nl) + aco2(nl)
       gtmp = (1.0 + 118.3*ag*rmp)**0.45
       gtmp2 = -1.41*ag*rmp
       garg = gtmp2/gtmp
c   Water Vapor
       wtmp = (1.0+20.07*awv(nl)*wvapor*rm)**0.45
       wtmp2 = -0.2385*awv(nl)*wvapor*rm
       warg = wtmp2/wtmp
       Tgas(nl) = exp(oarg+garg+warg)
      enddo
c
c  Compute clear sky transmittances
      call clrtrans(cosunz,rm,rmp,rmo,ws,relhum,am,Vi,
     *              ta,wa,asym,Tdclr,Tsclr)
      do nl = 1,nlt
       Foinc = Fobar(nl)*daycor*cosunz
c    Direct irradiance 
       Edir = Foinc*Tgas(nl)*Tdclr(nl)
c    Diffuse irradiance
       Edif = Foinc*Tgas(nl)*Tsclr(nl)
c    Spectral components
       Edclr(nl) = Edir
       Esclr(nl) = Edif
      enddo    !end clear nl loop
c  Compute cloudy transmittances
c      call slingo(rmu0,cldtc(ic,jc),rlwp(ic,jc),cdre(ic,jc),
c     *            Tdcld,Tscld)
      call slingo_verbose(rmu0,cldtau,clwp,re,
     *            Tdcld,Tscld)
      do nl = 1,nlt
       Foinc = Fobar(nl)*daycor*cosunz
c    Direct irradiance 
       Edir = Foinc*Tgas(nl)*Tdcld(ica(nl))
c    Diffuse irradiance
       Edif = Foinc*Tgas(nl)*Tscld(ica(nl))
c    Spectral components
       Edcld(nl) = Edir
       Escld(nl) = Edif
      enddo   !end cloudy nl loop
c
c  Sum clear and cloudy by percent
      ccov1 = cov*0.01  !convert from percent to fraction
      do nl = 1,nlt
       Ed(nl) = Edclr(nl)
       Es(nl) = Esclr(nl) 
       Ed(nl) = (1.0-ccov1)*Edclr(nl) + ccov1*Edcld(nl)
       Es(nl) = (1.0-ccov1)*Esclr(nl) + ccov1*Escld(nl)
      enddo
c
c   Total short-wave (W/m2)
       sirrsw = sirrsw + Edir + Edif
c
5     continue
      return
      end
