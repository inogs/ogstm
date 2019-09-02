      subroutine sfcsolz(iyr,iday,ihr,solz)
c
c  Computes solar zenith angle at surface given location 
c  and time.  
c
#include "define.h"
#include "rad.h"
#define RAD2HR 1
      save
      real solz(im,jm)
      data ifst /0/
c
c  Set up
      if (ifst .eq. 0)then
       call localvec
       nstps1 = 5
       delmin = 60.0/float(nstps1-1)
#if RAD2HR
c  compute irradiance every 2 hours
       nhr = 2
       nstps = nstps1*nhr - 1
       hrsec = 1.0/(float(nhr)*3600.0)
       hrsrt = 1.0
       hrend = 1.0
       iprt = 13
#else
       nstps = nstps1
       hrsec = 1.0/3600.0
       hrsrt = 1.0
       hrend = 0.0
       iprt = 12
#endif
       delh = delmin/60.0         !time interval in fraction of hour
       delx = delh*3600.0  !integration interval in seconds
       ifst = 1
      endif
c
c  Integrate to obtain mean cosine solar zenith angle
      solz = -1.0
      hr = float(ihr)

      do jl = 1,jm
        do il =1,im

           gmt = hr - hrsrt
           csza = 0.0

           do it = 1,nstps,(nstps-1)
             call sunmod(rad,iday,iyr,gmt,il,jl,sunz,rs)
             csza = csza + 0.5*cos(sunz/rad)
             gmt = hr + hrend
           enddo    !end of it loop

           gmt = hr - hrsrt

           do it = 2,nstps-1
             gmt = gmt + delh
             call sunmod(rad,iday,iyr,gmt,il,jl,sunz,rs)
             csza = csza + cos(sunz/rad)
           enddo   !end of it loop

           cosunz = csza*delx*hrsec*rad
           rsza = acos(cosunz/rad)
           sza = rsza*rad
           solz(il,jl) = min(sza,90.0)
           solz(il,jl) = max(solz(il,jl),0.0)

       enddo
      enddo
c
      return
      end
c
c **********************************************************************
      subroutine localvec
c
c  Create arrays of up, north, and east vectors for fixed locations
c  corresponding to these nwater indicators.  
c
#include "define.h"
#include "rad.h"
      real no
      common /bvec/ up(im,jm,3),no(im,jm,3),ea(im,jm,3)
c
c  Compute local east, north, and vertical vectors assuming longitude
c  is always 0
      do nv = 1,3
       do jl = 1,jm
         do il = 1,im
           no(il,jl,nv) = 0.0
           ea(il,jl,nv) = 0.0
           up(il,jl,nv) = 0.0
        enddo
       enddo
      enddo
c
c  Convert geodetic lat/lon to Earth-centered, earth-fixed (ECEF)
c  vector (geodetic unit vector)
      alon1h = alon0 - 0.5*dlon  !alat0 is bottom of 1st grid
      xlon1 = alon1h
      alat1h = alat0 - 0.5*dlat  !alat0 is bottom of 1st grid
      ylat1 = alat1h

      do jl = 1,jm
        do il = 1,im

         xlon = xlon1 + dlon * real(il,4)
         rlon = xlon/rad
         cosx = cos(rlon)
         sinx = sin(rlon)

         ylat = ylat1 + dlat * real(jl,4)
         rlat = ylat/rad
         cosy = cos(rlat)
         siny = sin(rlat)

         up(il,jl,1) = cosy*cosx
         up(il,jl,2) = cosy*sinx
         up(il,jl,3) = siny
c
c  Compute the local East and North unit vectors
         upxy = sqrt(up(il,jl,1)*up(il,jl,1)+up(il,jl,2)*up(il,jl,2))
         ea(il,jl,1) = -up(il,jl,2)/upxy
         ea(il,jl,2) = up(il,jl,1)/upxy
         no(il,jl,1) = up(il,jl,2)*ea(il,jl,3) - up(il,jl,3)*ea(il,jl,2)  !cross product
         no(il,jl,2) = up(il,jl,3)*ea(il,jl,1) - up(il,jl,1)*ea(il,jl,3)
         no(il,jl,3) = up(il,jl,1)*ea(il,jl,2) - up(il,jl,2)*ea(il,jl,1)
         enddo
      enddo

      write(6,*)'Local pointing vectors'
      il = 2
      jl = 2
      write(6,*)'no(2,2) = ',(no(il,jl,nv),nv=1,3)
      write(6,*)'ea(2,2) = ',(ea(il,jl,nv),nv=1,3)
      write(6,*)'up(2,2) = ',(up(il,jl,nv),nv=1,3)
c
      return
      end
