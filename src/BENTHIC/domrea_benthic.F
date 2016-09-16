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
      CHARACTER(LEN=11) B1Bmaskfile
      CHARACTER(LEN=50) filename
      CHARACTER(LEN=3), DIMENSION(7) :: var_nc
      CHARACTER(LEN=5) nomevar01
      LOGICAL B
      
! -------------------
! To read only one meshmask.nc

      maskfile        = 'meshmask.nc'
      call EXISTVAR(maskfile,'e3t_ben',B)
      if (B) then
         CALL readnc_slice_double_ben(maskfile,'e3t_ben', e3t_ben )
      else
         if(lwp) write(*,*) 'e3t_ben not found'
      endif      

      call EXISTVAR(maskfile,'e3w_ben',B)
      if (B) then
         CALL readnc_slice_double_ben(maskfile,'e3w_ben', e3w_ben )
      else
         if(lwp) write(*,*) 'e3w_ben not found'
      endif      

      call EXISTVAR(maskfile,'gdept_ben',B)
      if (B) then
         CALL readmask_double_1d_ben(maskfile,'gdept_ben', gdept_ben )
      else
         if(lwp) write(*,*) 'gdept_ben not found'
      endif      

      call EXISTVAR(maskfile,'gdepw_ben',B)
      if (B) then
         CALL readmask_double_1d_ben(maskfile,'gdepw_ben', gdepw_ben )
      else
         if(lwp) write(*,*) 'gdepw_ben not found'
      endif      


      call EXISTVAR(maskfile,'calcbenmask',B)
      if (B) then
      CALL readnc_slice_int_2d_ben (maskfile,'calcbenmask', calcbenmask )
      else
         if(lwp) write(*,*) 'calcbenmask not found'
      endif





!      B1Bmaskfile        = 'B1Bmask.nc'
!      call EXISTVAR(B1Bmaskfile,'calcbenflag',B)
!      if (B) then
!         CALL readmask_double_1d_ben(B1Bmaskfile,'calcbenflag',calcbenflag )
!      else
!         if(lwp) write(*,*) 'calcbenflag not found in ', B1Bmaskfile
!      endif      

      END SUBROUTINE domrea_benthic
