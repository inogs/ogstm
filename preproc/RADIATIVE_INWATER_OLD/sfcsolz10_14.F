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
