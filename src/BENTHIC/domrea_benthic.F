      SUBROUTINE domrea_benthic
!---------------------------------------------------------------------

! ======================
      USE calendar
      USE myalloc
      USE myalloc_mpp
      USE BIOBEN_mem
      USE TIME_MANAGER

      IMPLICIT NONE

! local declarations
! ==================

      INTEGER ii, jj, kk, jn, iinew, jjnew, iiend, jjend
      INTEGER ierr

      CHARACTER(LEN=11) maskfile
      CHARACTER(LEN=50) filename
      CHARACTER(LEN=3), DIMENSION(7) :: var_nc
      CHARACTER(LEN=5) nomevar01
      LOGICAL B
      
! -------------------
! To read only one meshmask.nc

      maskfile        = 'meshmask.nc'
      CALL readnc_slice_double (maskfile,'e3t_ben', e3t_ben )
      CALL readnc_slice_double (maskfile,'e3w_ben', e3w_ben )


      END SUBROUTINE domrea_benthic
