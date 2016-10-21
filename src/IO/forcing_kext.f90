      SUBROUTINE forcings_KEXT(datestring)
!---------------------------------------------------------------------


       USE myalloc
       USE myalloc_mpp
       USE OPT_mem
       USE TIME_MANAGER
       IMPLICIT NONE

      character(LEN=17), INTENT(IN) ::  datestring

! local declarations
! ==================
      REAL(8) sec,zweigh
      integer Before, After
      INTEGER iswap

      forcing_kext_partTime = MPI_WTIME()

      sec=datestring2sec(DATEstring)
      call TimeInterpolation(sec,TC_LEX, BEFORE, AFTER, zweigh)

      iswap  = 0


! ----------------------- INITIALIZATION -------------
      IF (datestring.eq.DATESTART) then

          CALL LOAD_KEXT(TC_LEX%TimeStrings(TC_LEX%Before)) ! CALL dynrea(iperm1)
          iswap = 1
          call swap_KEXT


        CALL LOAD_KEXT(TC_LEX%TimeStrings(TC_LEX%After)) !CALL dynrea(iper)


      ENDIF





! --------------------------------------------------------
! and now what we have to DO at every time step
! --------------------------------------------------------

! check the validity of the period in memory

      IF (BEFORE.ne.TC_LEX%Before) then
         TC_LEX%Before = BEFORE
         TC_LEX%After  = AFTER

         call swap_KEXT
         iswap = 1


          CALL LOAD_KEXT(TC_LEX%TimeStrings(TC_LEX%After))

          IF(lwp) WRITE (numout,*) ' Extinction factor DATA READ for Time = ', TC_LEX%TimeStrings(TC_LEX%After)
!      ******* LOADED NEW FRAME *************
      END IF





! compute the DATA at the given time step

      SELECT CASE (nsptint)
           CASE (0)  !  ------- no time interpolation
!      we have to initialize DATA IF we have changed the period
              IF (iswap.eq.1) THEN
                 zweigh = 1.0
                 call ACTUALIZE_KEXT(zweigh)! initialize now fields with the NEW DATA READ
              END IF

          CASE (1) ! ------------linear interpolation ---------------
             call ACTUALIZE_KEXT(zweigh)
      END SELECT


      forcing_kext_partTime = MPI_WTIME() - forcing_kext_partTime
      forcing_kext_TotTime  = forcing_kext_TotTime + forcing_kext_partTime



      RETURN
      END

! ******************************************************
!     SUBROUTINE LOAD_KEXT(datestring)
!     loads FORCINGS/KextF_yyyy0107-00:00:00.nc in OPT_mem.kextIO
! ******************************************************
      SUBROUTINE LOAD_KEXT(datestring)
! ======================
      USE calendar
      USE myalloc
      USE myalloc_mpp
      USE OPT_mem
      USE TIME_MANAGER
      USE BC_mem

      IMPLICIT NONE

      CHARACTER(LEN=17), INTENT(IN) :: datestring

      character(LEN=35) nomefile

      nomefile='FORCINGS/KextF_yyyy0107-00:00:00.nc'


!     Starting I/O
!    **********************************************************
      nomefile = 'FORCINGS/KextF_'//datestring//'.nc'

      if(lwp) write(*,'(A,I4,A,A)') "LOAD_KEXT --> I am ", rank, " starting reading forcing fields from ", nomefile

      call readnc_slice_float_2d(nomefile,'kextfact',buf2)
       kextIO(:,:,2) = buf2*tmask(:,:,1)
      

      END SUBROUTINE LOAD_KEXT





! ******************************************************
!     SUBROUTINE ACTUALIZE_PHYS(zweigh)
!     performs time interpolation
!     x(1)*(1-zweigh) + x(2)*zweigh
! ******************************************************
      SUBROUTINE ACTUALIZE_KEXT(zweigh)
         USE myalloc
         USE OPT_mem
         IMPLICIT NONE

         REAL(8), INTENT(IN) :: zweigh

        ! local
        INTEGER jj,ji

        DO ji=1,jpi
          DO jj=1,jpj
             kef(jj,ji) = ( (1. - zweigh) * kextIO(jj,ji,1)+ zweigh     * kextIO(jj,ji,2) )
          END DO
        END DO


      END SUBROUTINE ACTUALIZE_KEXT




! *************************************************************
!      SUBROUTINE SWAP_KEXT
! *    copia l'indice 2 nell'indice 1
! *************************************************************

      SUBROUTINE swap_KEXT
         USE myalloc
         USE OPT_mem
         IMPLICIT NONE


         INTEGER jj,ji

            DO ji=1,jpi
              DO jj=1,jpj
                kextIO(jj,ji,1) = kextIO(jj,ji,2)
              END DO
            END DO


      END SUBROUTINE swap_KEXT


