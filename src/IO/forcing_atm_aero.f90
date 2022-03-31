      SUBROUTINE forcings_atm_aero(datestring)
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
      call TimeInterpolation(sec,TC_OPTAERO, BEFORE, AFTER, zweigh)

      iswap  = 0


! ----------------------- INITIALIZATION -------------
      IF (datestring.eq.DATESTART) then

          CALL LOAD_atm_aero(TC_OPTAERO%TimeStrings(TC_OPTAERO%Before))
          iswap = 1
          call swap_atm_aero


        CALL LOAD_atm_aero(TC_OPTAERO%TimeStrings(TC_OPTAERO%After))


      ENDIF





! --------------------------------------------------------
! and now what we have to DO at every time step
! --------------------------------------------------------

! check the validity of the period in memory

      IF (BEFORE.ne.TC_OPTAERO%Before) then
         TC_OPTAERO%Before = BEFORE
         TC_OPTAERO%After  = AFTER

         call swap_atm_aero
         iswap = 1


          CALL LOAD_atm_aero(TC_OPTAERO%TimeStrings(TC_OPTAERO%After))

          IF(lwp) WRITE (numout,*) ' Extinction factor DATA READ for Time = ', TC_OPTAERO%TimeStrings(TC_OPTAERO%After)
!      ******* LOADED NEW FRAME *************
      END IF





! compute the DATA at the given time step

      SELECT CASE (nsptint)
           CASE (0)  !  ------- no time interpolation
!      we have to initialize DATA IF we have changed the period
              IF (iswap.eq.1) THEN
                 zweigh = 1.0
                 call ACTUALIZE_atm_aero(zweigh)! initialize now fields with the NEW DATA READ
              END IF

          CASE (1) ! ------------linear interpolation ---------------
             call ACTUALIZE_atm_aero(zweigh)
      END SELECT


      forcing_kext_partTime = MPI_WTIME() - forcing_kext_partTime
      forcing_kext_TotTime  = forcing_kext_TotTime + forcing_kext_partTime



      RETURN
      END

! ******************************************************
!     SUBROUTINE LOAD_atm_aero(datestring)
!     loads OPTICS/atm_yyyy0107-00:00:00.nc
! ******************************************************
      SUBROUTINE LOAD_atm_aero(datestring)
! ======================
      USE calendar
      USE myalloc
      USE OPT_mem
      USE TIME_MANAGER
      USE BC_mem

      IMPLICIT NONE

      CHARACTER(LEN=17), INTENT(IN) :: datestring
      integer jk
      character(LEN=32) nomefile

      nomefile='OPTICS/aero_yyyy0107-00:00:00.nc'


!     Starting I/O
!    **********************************************************
      nomefile = 'OPTICS/aero.'//datestring//'.nc'

      if(lwp) write(*,'(A,I4,A,A)') "LOAD_KEXT --> I am ", myrank, " starting reading atmospheric fields from ", nomefile

       call readnc_slice_float(nomefile,'taua',buf,0)
       do jk=1,33
         tauaIO(jk,:,:,2) = buf(jk,:,:)*tmask(1,:,:)
       enddo

       call readnc_slice_float(nomefile,'asymp',buf,0)
       do jk=1,33
          asympIO(jk,:,:,2) = buf(jk,:,:)*tmask(1,:,:)
       enddo

       call readnc_slice_float(nomefile,'ssalb',buf,0)
       do jk=1,33
          ssalbIO(jk,:,:,2) = buf(jk,:,:)*tmask(1,:,:)
       enddo




      

      END SUBROUTINE LOAD_atm_aero





! ******************************************************
!     SUBROUTINE ACTUALIZE_atm_aero(zweigh)
!     performs time interpolation
!     x(1)*(1-zweigh) + x(2)*zweigh
! ******************************************************
      SUBROUTINE ACTUALIZE_atm_aero(zweigh)
         USE myalloc
         USE OPT_mem
         IMPLICIT NONE

         double precision, INTENT(IN) :: zweigh
         double precision Umzweigh

         Umzweigh =1.0 - zweigh

        taua  = Umzweigh*  tauaIO(:,:,:,1) + zweigh*  tauaIO(:,:,:,2)
        asymp = Umzweigh* asympIO(:,:,:,1) + zweigh* asympIO(:,:,:,2)
        ssalb = Umzweigh* ssalbIO(:,:,:,1) + zweigh* ssalbIO(:,:,:,2)





      END SUBROUTINE ACTUALIZE_atm_aero




! *************************************************************
!      SUBROUTINE SWAP_KEXT
! *    copia l'indice 2 nell'indice 1
! *************************************************************

      SUBROUTINE swap_atm_aero
         USE myalloc
         USE OPT_mem
         IMPLICIT NONE

         tauaIO(:,:,:,1)  = tauaIO(:,:,:,2)
        asympIO(:,:,:,1) = asympIO(:,:,:,2)
        ssalbIO(:,:,:,1) = ssalbIO(:,:,:,2)






      END SUBROUTINE swap_atm_aero


