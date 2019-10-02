      subroutine sfcsolz10_14(iyr,iday,ihr,solz)
      USE myalloc
      USE OPT_mem
      IMPLICIT NONE
!  Computes solar zenith angle at surface given location 
!  and time.  
!
      integer, intent(IN)   :: iyr, iday, ihr
      double precision, intent(OUT)  :: solz(jpj,jpi)

! local variables
      integer               :: ji, jj, it
      double precision      :: hr, gmt, csza, sunz, rs
      double precision      :: cosunz, rsza, sza
!
!  Set up
      if (ifst) then
       call localvec
!      nstps1 = 5
       delmin = 60.0d0/real(nstps1-1,8)

       nstps = nstps1
       hrsec = 1.0d0/3600.0d0
!      original values integration is done the hors before
!      hrsrt = 1.0d0
!      hrend = 0.0d0
!      we integrate for the following our
       hrsrt = 0.0d0
       hrend = 4.0d0

       delh = delmin/60.0d0         !time interval in fraction of hour
       delx = delh*3600.0d0  !integration interval in seconds
       ifst = .FALSE.
      endif
!
!  Integrate to obtain mean cosine solar zenith angle
      solz = -1.0d0
      hr = real(ihr,8)

      do ji = 1,jpi
        do jj =1,jpj

           gmt = hr - hrsrt
           csza = 0.0d0

           do it = 1,nstps,(nstps-1)
             call sunmod(rad,iday,iyr,gmt,jj,ji,sunz,rs)
             csza = csza + 0.5d0*cos(sunz/rad)
             gmt = hr + hrend
           enddo    !end of it loop

           gmt = hr - hrsrt

           do it = 2,nstps-1
             gmt = gmt + delh
             call sunmod(rad,iday,iyr,gmt,jj,ji,sunz,rs)
             csza = csza + cos(sunz/rad)
           enddo   !end of it loop

!          csza*delx is the integral, hrsec is the normalization factor 

           cosunz = csza*delx*hrsec*rad

           rsza = acos(cosunz/rad) ! average angle in radians
           sza = rsza*rad ! conversion in degrees

           solz(jj,ji) = min(sza,90.0d0)
           solz(jj,ji) = max(solz(jj,ji),0.0d0)

       enddo
      enddo
c
      return
      end
c
c **********************************************************************
      subroutine localvec
      USE myalloc, ONLY: jpi,jpj,glamt,gphit
      USE OPT_mem, ONLY: rad,no,ea,up

      IMPLICIT NONE

C local variables
      integer               :: ji, jj, nv
      double precision      :: xlon, rlon, cosx, sinx 
      double precision      :: ylat, rlat, cosy, siny
      double precision      :: upxy
c
c  Create arrays of up, north, and east vectors for fixed locations
c  corresponding to these nwater indicators.  
c
c
c  Compute local east, north, and vertical vectors assuming longitude
c  is always 0
      do nv = 1,3
       do ji = 1,jpi
        do jj = 1,jpj
           no(jj,ji,nv) = 0.0d0
           ea(jj,ji,nv) = 0.0d0
           up(jj,ji,nv) = 0.0d0
        enddo
       enddo
      enddo
c
c  Convert geodetic lat/lon to Earth-centered, earth-fixed (ECEF)
c  vector (geodetic unit vector)
c     alon1h = alon0 - 0.5*dlon  !alat0 is bottom of 1st grid
c     xlon1 = alon1h
c     alat1h = alat0 - 0.5*dlat  !alat0 is bottom of 1st grid
c     ylat1 = alat1h

      do ji = 1,jpi
        do jj = 1,jpj

         xlon = 0.0d0
C        xlon = glamt(jj,ji) ! commenting this line is equivalnet to  GMT_ALL
         rlon = xlon/rad
         cosx = cos(rlon)
         sinx = sin(rlon)

         ylat = gphit(jj,ji)
         rlat = ylat/rad
         cosy = cos(rlat)
         siny = sin(rlat)

         up(jj,ji,1) = cosy*cosx
         up(jj,ji,2) = cosy*sinx
         up(jj,ji,3) = siny
c
c  Compute the local East and North unit vectors
         upxy = sqrt(up(jj,ji,1)*up(jj,ji,1)+up(jj,ji,2)*up(jj,ji,2))
         ea(jj,ji,1) = -up(jj,ji,2)/upxy
         ea(jj,ji,2) = up(jj,ji,1)/upxy
         no(jj,ji,1) = up(jj,ji,2)*ea(jj,ji,3) - up(jj,ji,3)*ea(jj,ji,2)  !cross product
         no(jj,ji,2) = up(jj,ji,3)*ea(jj,ji,1) - up(jj,ji,1)*ea(jj,ji,3)
         no(jj,ji,3) = up(jj,ji,1)*ea(jj,ji,2) - up(jj,ji,2)*ea(jj,ji,1)
         enddo
      enddo

c     write(6,*)'Local pointing vectors'
c     ji = 2
c     jj = 2
c     write(6,*)'no(2,2) = ',(no(jj,ji,nv),nv=1,3)
c     write(6,*)'ea(2,2) = ',(ea(jj,ji,nv),nv=1,3)
c     write(6,*)'up(2,2) = ',(up(jj,ji,nv),nv=1,3)
c
      return
      end
