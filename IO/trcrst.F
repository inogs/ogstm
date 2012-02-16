      SUBROUTINE trcrst
!---------------------------------------------------------------------
!
!                       ROUTINE trcrst
!                     ******************
!
!  PURPOSE :
!  ---------
!     READ files for restart for passive tracer
!
!----------------------------------------------------------------------


       USE calendar
       USE myalloc
       USE myalloc_mpp
       USE TIME_MANAGER

       IMPLICIT NONE

!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER ii, jj, kk, jn, iinew, jjnew, iista, iiend, jjsta, jjend
      CHARACTER*29 filename
      CHARACTER*7  varname


#include "stafun.h"

! 0. initialisations
       filename = 'RST.20111231-15:30:00.N1p.nc'
       filename(29:29)=CHAR(0)
       varname = 'TRNN1p'
       varname(7:7) = CHAR(0)

!      qui leggeva il date.step

!     - To read only one meshmask.nc
      iista = nimpp
      iiend = MIN(nimpp + jpi - 1,jpiglo)
      jjsta = njmpp
      jjend = MIN(njmpp + jpj - 1,jpjglo)
      totglamt = 0.0 ! a cosa serve ?
      totgphit = 0.0
      rsttrn   = 0.0
      rsttrb   = 0.0

      IF(lwp) THEN
          WRITE(numout,*) ' '
          WRITE(numout,*) ' *** trcrst beginning of restart for'
          WRITE(numout,*) ' passive tracer'
          WRITE(numout,*) '   with the time TimeStepStart : ',TimeStepStart
          WRITE(numout,*) ' '
      ENDIF

      DO jn=1, jptra  ! global loop on tracers to read restart


         varname( 1:6)  = 'TRN'//ctrcnm(jn)(1:3)
         filename(1:28) = 'RST.'//DateStart//'.'//ctrcnm(jn)(1:3)//'.nc'

         CALL IOOGSNC(filename, varname, jpiglo, jpjglo, jpk, rsttrn(1,1,1))

         varname(1:6) = 'TRB'//ctrcnm(jn)(1:3)
         CALL IOOGSNC(filename, varname, jpiglo, jpjglo, jpk, rsttrb(1,1,1))


          do kk=1, jpk
           do jj=jjsta, jjend
             do ii=iista, iiend
               iinew = ii - nimpp + 1
               jjnew = jj - njmpp + 1
               trn(iinew, jjnew, kk, jn) = rsttrn(ii, jj, kk)
               trb(iinew, jjnew, kk, jn) = rsttrb(ii, jj, kk)
             enddo
           enddo
         enddo


! ********************   we put initial undef to 0
          do kk=1, jpk
           do jj=1, jpj
             do ii=1, jpi
               trn(ii, jj, kk, jn) = trn(ii, jj, kk, jn)*tmask(ii,jj,kk)
               trb(ii, jj, kk, jn) = trb(ii, jj, kk, jn)*tmask(ii,jj,kk)
            enddo
           enddo
         enddo

       END DO




      END SUBROUTINE trcrst

