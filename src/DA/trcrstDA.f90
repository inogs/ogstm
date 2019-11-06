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
      INTEGER jn, jDA, shift, AssimilationLevels
      CHARACTER(LEN=45) filename
      CHARACTER(LEN=6)  varname
      CHARACTER(LEN=2)  HOUR
      CHARACTER(LEN=5)  OBStype


      shift = 0

      HOUR = datestring(10:11)
      if(HOUR.eq.'12') then
        AssimilationLevels = AssimilationLevels_sat
      elseif(HOUR.eq.'13') then
        AssimilationLevels = AssimilationLevels_float
      endif

      DO jn=1, jptra  ! global loop on tracers to read restart


      if (.not.isaDAvar(ctrcnm(jn))) CYCLE


         varname  = 'TRN'//ctrcnm(jn)
         filename = 'DA__FREQ_1/RST_after.'//datestring//'.'//trim(ctrcnm(jn))//'.nc'
         if (lwp) write(*,*) 'reading ', filename
         CALL readnc_slice_floatDA(filename,varname,AssimilationLevels,trn(:,:,:,jn), shift)



! ********************   we put initial undef to 0
          trb(:,:,:,jn) = trn(:,:,:,jn) * tmask; !! ACHTUNG !!!
          trn(:,:,:,jn) = trn(:,:,:,jn) * tmask;
      ENDDO


      END SUBROUTINE trcrstDA

      SUBROUTINE readnc_slice_floatDA(fileNetCDF,varname,AssimilationLevels, M, shift)
      USE myalloc
      USE netcdf
      ! USE DA_mem, ONLY : AssimilationLevels_sat,AssimilationLevels_float
      implicit none


      character,intent(in) :: fileNetCDF*(*) ,varname*(*)
      integer, intent(in)  :: shift, AssimilationLevels
      double precision,intent(inout) ::  M(jpk,jpj,jpi)
      
      real,allocatable,dimension(:,:,:) :: copy_in
      integer ncid, stat, VARid,i,j,k
      integer counter
      integer thecount(4), start(4)

      
      allocate(copy_in(jpi,jpj,AssimilationLevels))
      counter = 0
      start    = (/nimpp+shift, njmpp,  1,  1/)
      thecount = (/jpi,           jpj, AssimilationLevels, 1/)

      stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)  
       call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_inq_varid (ncid, varname, VARid)
       call handle_err2(stat, fileNetCDF,varname)
       call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_get_var (ncid,VARid,copy_in,start, thecount)
      
      call handle_err2(stat, fileNetCDF,varname)        
      call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_close(ncid)                           
      call handle_err1(stat, counter,FileNetCDF)

      DO i=1,jpi
        DO j=1,jpj
          DO k=1,AssimilationLevels
            M(k,j,i) = real(copy_in(i,j,k),8)
          ENDDO
        ENDDO
      ENDDO   
      
      deallocate(copy_in)

      END SUBROUTINE readnc_slice_floatDA

