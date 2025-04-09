      SUBROUTINE forcings_atm_clim(datestring)
!---------------------------------------------------------------------


       USE myalloc
       ! epascolo USE myalloc_mpp
       USE OPT_mem
       USE TIME_MANAGER
       USE mpi
       IMPLICIT NONE

      character(LEN=17), INTENT(IN) ::  datestring

! local declarations
! ==================
      double precision sec,zweigh
      integer Before, After
      INTEGER iswap

      forcing_kext_partTime = MPI_WTIME()

      sec=datestring2sec(DATEstring)
      call TimeInterpolation(sec,TC_OPTCLIM, BEFORE, AFTER, zweigh)

      iswap  = 0


! ----------------------- INITIALIZATION -------------
      IF (datestring.eq.DATESTART) then

          CALL LOAD_climatm(TC_OPTCLIM%TimeStrings(TC_OPTCLIM%Before))
          iswap = 1
          call swap_climatm


        CALL LOAD_climatm(TC_OPTCLIM%TimeStrings(TC_OPTCLIM%After))


      ENDIF





! --------------------------------------------------------
! and now what we have to DO at every time step
! --------------------------------------------------------

! check the validity of the period in memory

      IF (BEFORE.ne.TC_OPTCLIM%Before) then
         TC_OPTCLIM%Before = BEFORE
         TC_OPTCLIM%After  = AFTER

         call swap_KEXT
         iswap = 1


          CALL LOAD_climatm(TC_OPTCLIM%TimeStrings(TC_OPTCLIM%After))

          IF(lwp) WRITE (numout,*) ' Extinction factor DATA READ for Time = ', TC_OPTCLIM%TimeStrings(TC_OPTCLIM%After)
!      ******* LOADED NEW FRAME *************
      END IF





! compute the DATA at the given time step

      SELECT CASE (nsptint)
           CASE (0)  !  ------- no time interpolation
!      we have to initialize DATA IF we have changed the period
              IF (iswap.eq.1) THEN
                 zweigh = 1.0
                 call ACTUALIZE_climatm(zweigh)! initialize now fields with the NEW DATA READ
              END IF

          CASE (1) ! ------------linear interpolation ---------------
             call ACTUALIZE_climatm(zweigh)
      END SELECT


      forcing_kext_partTime = MPI_WTIME() - forcing_kext_partTime
      forcing_kext_TotTime  = forcing_kext_TotTime + forcing_kext_partTime



      RETURN
      END

! ******************************************************
!     SUBROUTINE LOAD_climatm(datestring)
!     loads OPTICS/atm_yyyy0107-00:00:00.nc
! ******************************************************
      SUBROUTINE LOAD_climatm(datestring)
! ======================
      USE calendar
      USE myalloc
      USE OPT_mem
      USE TIME_MANAGER
      USE BC_mem

      IMPLICIT NONE

      CHARACTER(LEN=17), INTENT(IN) :: datestring

      character(LEN=35) nomefile

      nomefile='OPTICS/climatm_yyyy0107-00:00:00.nc'


!     Starting I/O
!    **********************************************************
      nomefile = 'OPTICS/climatm.'//datestring//'.nc'

      if(lwp) write(*,'(A,I4,A,A)') "LOAD_KEXT --> I am ", myrank, " starting reading atmospheric fields from ", nomefile


       call readnc_slice_float_2d(nomefile,'cdrem',buf2,0)
       cdremIO(:,:,2) = buf2*tmask(1,:,:)

       call readnc_slice_float_2d(nomefile,'cldtcm',buf2,0)
       cldtcmIO(:,:,2) = buf2*tmask(1,:,:)


      

      END SUBROUTINE LOAD_climatm





! ******************************************************
!     SUBROUTINE ACTUALIZE_climatm(zweigh)
!     performs time interpolation
!     x(1)*(1-zweigh) + x(2)*zweigh
! ******************************************************
      SUBROUTINE ACTUALIZE_climatm(zweigh)
         USE myalloc
         USE OPT_mem
         IMPLICIT NONE

         double precision, INTENT(IN) :: zweigh


        call actualize(zweigh,cdremIO,cdrem)
        call actualize(zweigh,cldtcmIO,cldtcm)





      END SUBROUTINE ACTUALIZE_climatm




! *************************************************************
!      SUBROUTINE SWAP_KEXT
! *    copia l'indice 2 nell'indice 1
! *************************************************************

      SUBROUTINE swap_climatm
         USE myalloc
         USE OPT_mem
         IMPLICIT NONE

        call swap(cdremIO)
        call swap(cldtcmIO)




      END SUBROUTINE swap_climatm


