      subroutine slingo(rmu0,cldtau,clwp,cre,Tcd,Tcs)
c
c  Slingo's (1989) Delta-Eddington approximation for the two-
c  stream equations applied to clouds.  
c  Inputs:
c       rmu0    Kasten's approx for cosine of solar zenith angle
c       cldtau  cloud optical thickness (at 0.6 um)
c       clwp    liquid water path in cloud (g/m2)
c       cre     cloud droplet effective radius (um)
c  Outputs
c       ica     index for translating cloud arrays to light arrays
c       Tcd     spectral transmittance for downwelling direct irradiance
c       Tcs     spectral transmittance for downwelling diffuse irradiance
c  Internal
c       Tcu     spectral transmittance for upwelling irradiance
c
c
#include "comlte.h"
      save ifst,U1
      parameter(ncld=24)
      common /bslg1/ asl(ncld),bsl(ncld),csl(ncld),dsl(ncld)
      common /bslg2/ esl(ncld),fsl(ncld)
      common /bica/ ica(nlt)
      real rnl1(ncld),rnl2(ncld)
      real Tcd(ncld),Tcs(ncld)
      real Tcu(ncld)
      data ifst /0/
c
c  Get one-time-only constants
      if (ifst .eq. 0)then
       U1 = 7.0/4.0
       call rdslingo(rnl1,rnl2)    !reads cloud parameter file
c   Indices to relate cloud parameters to other light parameters
       do nl = 1,nlt
        do nc = 1,ncld
         lamcld = nint(rnl2(nc)*1000.0)
         if (lam(nl) .lt. lamcld)then
          ica(nl) = nc
          go to 5
         endif
        enddo
5       continue
       enddo
c       ica(nlt) = ncld
       do nc = 1,ncld
        Tcd(nc) = 0.0
        Tcs(nc) = 0.0
        Tcu(nc) = 0.0
       enddo
       ifst = 1
      endif
c      
c  Compute re as funtion of cldtau and LWP according to eq. 1 in 
c  Slingo.
c   tau is derived at this wavelength (0.6 um) in the ISCCP data set
c      re = clwp*bsl(9)/(cldtau - clwp*asl(9))
c      re = min(re,15.0)  !block high re -- produces excessive direct
c  Changes to the ISCCP-D2 data set make this relationship untenable
c  (excessive re's are derived).  Instead choose a fixed re of 10 um
c  for ocean (Kiehl et al., 1998 -- J. Clim.)
c       re = 10.0
c  Paper by Han et al., 1994 (J.Clim.) show mean ocean cloud radius
c  = 11.8 um
c       re = 11.8
c  Mean of Kiehl and Han
      re = (10.0+11.8)/2.0
c
c  Compute spectral cloud characteristics
c   If MODIS re is available use it; otherwise use parameterized re above
      if (cre .ge. 0.0)then   !use MODIS re
       re = cre
      endif
      izero = 0
      do nc = 1,ncld
       tauc = clwp*(asl(nc)+bsl(nc)/re)
       oneomega = csl(nc) + dsl(nc)*re
       omega = 1.0 - oneomega
       g = esl(nc) + fsl(nc)*re
       b0 = 3.0/7.0*(1.0-g)
       bmu0 = 0.5 - 0.75*rmu0*g/(1.0+g)
       f = g*g
       U2 = U1*(1.0-((1.0-omega)/(7.0*omega*b0)))
       U2 = max(U2,0.0)
       alpha1 = U1*(1.0-omega*(1.0-b0))
       alpha2 = U2*omega*b0
       alpha3 = (1.0-f)*omega*bmu0
       alpha4 = (1.0-f)*omega*(1.0-bmu0)
       sqarg = alpha1*alpha1 - alpha2*alpha2
       sqarg = max(sqarg,1.0E-9)
C bioptimod  sqarg = max(sqarg,1.0E-17)
       eps = sqrt(sqarg)
       rM = alpha2/(alpha1+eps)
       E = exp(-eps*tauc)
       val1 = 1.0 - omega*f
       val2 = eps*eps*rmu0*rmu0
       rnum = val1*alpha3 - rmu0*(alpha1*alpha3+alpha2*alpha4)
       rden = val1*val1 - val2
c bioptimod added  + 0.00001 to avoid NaN
       gama1 = rnum/(rden + 0.00001)
       rnum = -val1*alpha4 - rmu0*(alpha1*alpha4+alpha2*alpha3)
c bioptimod added  + 0.00001 to avoid NaN
       gama2 = rnum/(rden + 0.00001)
       Tdb = exp(-val1*tauc/rmu0)
       val3 = 1.0 - E*E*rM*rM
c bioptimod added  + 0.00001 to avoid NaN
       Rdif = rM*(1.0-E*E)/(val3 + 0.00001)
c bioptimod added  + 0.00001 to avoid NaN
       Tdif = E*(1.0-rM*rM)/(val3 + 0.00001)
c       Rdir = -gama2*Rdif - gama1*Tdb*Tdif + gama1
       Tdir = -gama2*Tdif - gama1*Tdb*Rdif + gama2*Tdb
c       Tdir = max(Tdir,0.0)
       Tcd(nc) = Tdb
       Tcs(nc) = Tdir
       if (Tcs(nc) .lt. 0.0)izero=1
c       Tcu(nc) = Tdif
      enddo
      if (izero .eq. 1)then    !block negative diffuse irrads
       do nc = 1,ncld
        Tcs(nc) = 0.0
       enddo
      endif
c
      return
      end
c
c **********************************************************************
      subroutine rdslingo(rnl1,rnl2)
c
c  Reads cloud parameters by Slingo (1989).
c
      parameter(ncld=24)
      character*50 title
      real*4 rn1,rn2,a4,b4,c4,d4,e4,f4
      real rnl1(ncld),rnl2(ncld)
      common /bslg1/ asl(ncld),bsl(ncld),csl(ncld),dsl(ncld)
      common /bslg2/ esl(ncld),fsl(ncld)
c
      open(15,file='bcs/slingo.dat',status='old',form='formatted')
      do n = 1,3
       read(15,'(a50)')title
      enddo
      do n = 1,ncld
       read(15,10)rn1,rn2,a4,b4,e4,f4,c4,d4
       rnl1(n) = rn1
       rnl2(n) = rn2
       asl(n) = a4*0.01
       bsl(n) = b4
       csl(n) = c4
       dsl(n) = d4
       esl(n) = e4
       fsl(n) = f4*0.001
      enddo
      close(15)
c
10    format(2f5.2,3x,2f6.3,2f6.3,2f12.8)
c     format(2f5.2,3x,2f6.3,2f6.3,1pe9.2,1pe8.2)
      return
      end
