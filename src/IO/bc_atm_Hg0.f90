MODULE mod_Hgbc

USE myalloc
USE BIO_mem
USE TIME_MANAGER
USE calendar
USE mpi

implicit NONE

contains

      SUBROUTINE BC_atm_Hg0(datestring)
       
       IMPLICIT NONE

      character(LEN=17), INTENT(IN) ::  datestring

! local declarations
! ==================
      double precision sec,zweigh
      integer Before, After
      INTEGER iswap

      bc_atm_Hg0_partTime = MPI_WTIME()

      sec=datestring2sec(DATEstring)
      call TimeInterpolation(sec,TC_Hg0, BEFORE, AFTER, zweigh)

      iswap  = 0

! ----------------------- INITIALIZATION -------------
      IF (datestring.eq.DATESTART) then
!         scrivi qualche log sul numout

          CALL LOAD_atm_Hg0(TC_Hg0%TimeStrings(TC_Hg0%Before)) ! CALL trcrea(iperm1,bc)
          iswap = 1
          call swap_Hg0


        CALL LOAD_atm_Hg0(TC_Hg0%TimeStrings(TC_Hg0%After))    ! CALL trcrea(iper,bc)


      ENDIF

! --------------------------------------------------------
! and now what we have to DO at every time step
! --------------------------------------------------------

! check the validity of the period in memory
      IF (BEFORE.ne.TC_Hg0%Before) then
         TC_Hg0%Before = BEFORE
         TC_Hg0%After  = AFTER

         call swap_Hg0
         iswap = 1


          CALL LOAD_atm_Hg0(TC_Hg0%TimeStrings(TC_Hg0%After))

          IF(lwp) WRITE (numout,*) 'Hg0 atm  DATA READ for Time = ', TC_Hg0%TimeStrings(TC_Hg0%After)
!      ******* LOADED NEW FRAME *************
      END IF

! compute the DATA at the given time step

      SELECT CASE (nsptint)
           CASE (0)  !  ------- no time interpolation
!      we have to initialize DATA IF we have changed the period
              IF (iswap.eq.1) THEN
                 zweigh = 1.0
                 call ACTUALIZE_Hg0(zweigh)! initialize now fields with the NEW DATA READ
              END IF

          CASE (1) ! ------------linear interpolation ---------------
             call ACTUALIZE_Hg0(zweigh)

      END SELECT


      bc_atm_Hg0_partTime = MPI_WTIME()    - bc_atm_Hg0_partTime
      bc_atm_Hg0_TotTime = bc_atm_Hg0_TotTime + bc_atm_Hg0_partTime

      END SUBROUTINE BC_atm_Hg0




! ******************************************************
!     SUBROUTINE LOAD_atm_Hg0(datestring)
!     loads BC/atm_Hg0_yyyy0107-00:00:00.nc in BC_mem.Hg0_dtatrc(:,2,:)
! ******************************************************
      SUBROUTINE LOAD_atm_Hg0(datestring)
         
         

          IMPLICIT NONE

          CHARACTER(LEN=17), INTENT(IN) :: datestring

!         local
          character(LEN=31) nomefile

          nomefile='BC/atm_Hg0_yyyy0107-00:00:00.nc'


    !     Starting I/O
    !    **********************************************************
         nomefile(1:31)='BC/atm_Hg0_'//datestring//'.nc'
         if(lwp) write(*,'(A,I4,A,A)') "LOAD_atm_Hg0 --> I am ", myrank, " starting reading forcing fields from ", nomefile(1:31)

         call readnc_slice_float_2d(nomefile,'hg0atm',buf2,0)
         Hg0_IO(:,:,2) = buf2*tmask(1,:,:)
      


      END SUBROUTINE LOAD_atm_Hg0


! ******************************************************
!     SUBROUTINE ACTUALIZE_Hg0(zweigh)
!     performs time interpolation
!     x(1)*(1-zweigh) + x(2)*zweigh
! ******************************************************
      SUBROUTINE ACTUALIZE_Hg0(zweigh)
        
         IMPLICIT NONE

         double precision, INTENT(IN) :: zweigh

        ! local
        INTEGER jj,ji

        DO ji=1,jpi
          DO jj=1,jpj
             atm_Hg0(jj,ji) = ( (1. - zweigh) * Hg0_IO(jj,ji,1)+ zweigh     * Hg0_IO(jj,ji,2) )
          END DO
        END DO


      END SUBROUTINE ACTUALIZE_Hg0


! ****************************************************


      SUBROUTINE swap_Hg0
         
         IMPLICIT NONE


         INTEGER jj,ji

            DO ji=1,jpi
              DO jj=1,jpj
                Hg0_IO(jj,ji,1) = Hg0_IO(jj,ji,2)
              END DO
            END DO


      END SUBROUTINE swap_Hg0

END MODULE mod_Hgbc
