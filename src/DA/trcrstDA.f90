      SUBROUTINE trcrstDA(datestring)
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
       ! epascolo USE myalloc_mpp
       USE TIME_MANAGER
       USE DA_mem

       IMPLICIT NONE
       CHARACTER(LEN=17), INTENT(IN) :: datestring

!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER jn, jDA, shift
      CHARACTER(LEN=37) filename
      CHARACTER(LEN=6)  varname
      CHARACTER(LEN=43) bkpname


      shift = 0

      DO jn=1, jptra  ! global loop on tracers to read restart


      if (.not.isaDAvar(ctrcnm(jn))) CYCLE


         varname  = 'TRN'//ctrcnm(jn)
         filename = 'RESTARTS/RST.'//datestring//'.'//trim(ctrcnm(jn))//'.nc'
         if (lwp) write(*,*) 'reading ', filename
         CALL readnc_slice_float(filename,varname, trn(:,:,:,jn), shift)



! ********************   we put initial undef to 0
          trb(:,:,:,jn) = trn(:,:,:,jn) * tmask; !! ACHTUNG !!!
          trn(:,:,:,jn) = trn(:,:,:,jn) * tmask;
      ENDDO


      END SUBROUTINE trcrstDA

