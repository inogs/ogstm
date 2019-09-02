      subroutine slingo_verbose(rmu0,cldtau,clwp,cre,Tcd,Tcs)
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
       call rdslingo_verbose(rnl1,rnl2)    !reads cloud parameter file
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
      write(666,*) "re=", re 
      izero = 0
      do nc = 1,ncld
       tauc = clwp*(asl(nc)+bsl(nc)/re)
       if (nc .EQ. ica(1)) write(666,*) "clwp=", clwp
       if (nc .EQ. ica(1)) write(666,*) "tauc=", tauc
       if (nc .EQ. ica(1)) write(666,*) "asl(n)=", asl(nc)
       if (nc .EQ. ica(1)) write(666,*) "bsl(n)=", bsl(nc)
       oneomega = csl(nc) + dsl(nc)*re
       if (nc .EQ. ica(1)) write(666,*) "csl(nc)=", csl(nc)
       if (nc .EQ. ica(1)) write(666,*) "dsl(nc)=", dsl(nc)
       if (nc .EQ. ica(1)) write(666,*) "oneomega=", oneomega
       omega = 1.0 - oneomega
       g = esl(nc) + fsl(nc)*re
       if (nc .EQ. ica(1)) write(666,*)  "g=", g
       b0 = 3.0/7.0*(1.0-g)
       if (nc .EQ. ica(1)) write(666,*)  "b0=",b0
       bmu0 = 0.5 - 0.75*rmu0*g/(1.0+g)
       if (nc .EQ. ica(1)) write(666,*)  "bmu0=",bmu0
       f = g*g
       if (nc .EQ. ica(1)) write(666,*)  "f=",f
       U2 = U1*(1.0-((1.0-omega)/(7.0*omega*b0)))
       if (nc .EQ. ica(1)) write(666,*)  "U2=",U2
       U2 = max(U2,0.0)
       alpha1 = U1*(1.0-omega*(1.0-b0))
       if (nc .EQ. ica(1)) write(666,*)  "aplha1=",alpha1
       alpha2 = U2*omega*b0
       if (nc .EQ. ica(1)) write(666,*)  "aplha2=",alpha2
       alpha3 = (1.0-f)*omega*bmu0
       if (nc .EQ. ica(1)) write(666,*)  "aplha3=",alpha3
       alpha4 = (1.0-f)*omega*(1.0-bmu0)
       if (nc .EQ. ica(1)) write(666,*)  "aplha4=",alpha4
       sqarg = alpha1*alpha1 - alpha2*alpha2
       if (nc .EQ. ica(1)) write(666,*)  "sqarg=",sqarg
       sqarg = max(sqarg,1.0E-9)
c      sqarg = max(sqarg,1.0E-17)
       if (nc .EQ. ica(1)) write(666,*)  "max_sqarg=",sqarg
       eps = sqrt(sqarg)
       if (nc .EQ. ica(1)) write(666,*)  "eps=",eps
       rM = alpha2/(alpha1+eps)
       if (nc .EQ. ica(1)) write(666,*)  "rM=",rM
       E = exp(-eps*tauc)
       if (nc .EQ. ica(1)) write(666,*)  "E=",E
       val1 = 1.0 - omega*f
       if (nc .EQ. ica(1)) write(666,*)  "val1=",val1
       val2 = eps*eps*rmu0*rmu0
       if (nc .EQ. ica(1)) write(666,*)  "val2=",val2
       rnum = val1*alpha3 - rmu0*(alpha1*alpha3+alpha2*alpha4)
       if (nc .EQ. ica(1)) write(666,*)  "rnum=",rnum
       rden = val1*val1 - val2
       if (nc .EQ. ica(1)) write(666,*)  "rden=",rden
       gama1 = rnum/rden
       if (nc .EQ. ica(1)) write(666,*)  "gama1=",gama1
       rnum = -val1*alpha4 - rmu0*(alpha1*alpha4+alpha2*alpha3)
       if (nc .EQ. ica(1)) write(666,*)  "rnum=",rnum
       gama2 = rnum/rden
       if (nc .EQ. ica(1)) write(666,*)  "gama2=",gama2
       Tdb = exp(-val1*tauc/rmu0)
       if (nc .EQ. ica(1)) write(666,*)  "Tdb=",Tdb
       val3 = 1.0 - E*E*rM*rM
       if (nc .EQ. ica(1)) write(666,*)  "val3=",val3
       Rdif = rM*(1.0-E*E)/val3
       if (nc .EQ. ica(1)) write(666,*)  "Rdif=",Rdif
       Tdif = E*(1.0-rM*rM)/val3
       if (nc .EQ. ica(1)) write(666,*)  "Tdif=",Tdif
c       Rdir = -gama2*Rdif - gama1*Tdb*Tdif + gama1
       Tdir = -gama2*Tdif - gama1*Tdb*Rdif + gama2*Tdb
       if (nc .EQ. ica(1)) write(666,*)  "Tdir=",Tdir
c       Tdir = max(Tdir,0.0)
       Tcd(nc) = Tdb
       if (nc .EQ. ica(1)) write(666,*)  "Tcd=",Tcd(nc)
       Tcs(nc) = Tdir
       if (nc .EQ. ica(1)) write(666,*)  "Tcs=",Tcs(nc)
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
      subroutine rdslingo_verbose(rnl1,rnl2)
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
       write(666,*) "rnl1(n)=", rnl1(n)
       write(666,*) "rnl2(n)=", rnl2(n)
       write(666,*) "asl(n)=", asl(n)
       write(666,*) "bsl(n)=", bsl(n)
       write(666,*) "csl(n)=", csl(n)
       write(666,*) "dsl(n)=", dsl(n)
       write(666,*) "esl(n)=", esl(n)
       write(666,*) "fsl(n)=", fsl(n)
      enddo
      close(15)
c
10    format(2f5.2,3x,2f6.3,2f6.3,2f12.8)
      return
      end
