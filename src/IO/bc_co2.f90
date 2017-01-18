MODULE mod_cbc

USE myalloc
USE BIO_mem
USE TIME_MANAGER
USE calendar
USE mpi

implicit NONE

contains

      SUBROUTINE BC_CO2(datestring)
       
       IMPLICIT NONE

      character(LEN=17), INTENT(IN) ::  datestring

! local declarations
! ==================
      double precision sec,zweigh
      integer Before, After
      INTEGER iswap

      bc_co2_partTime = MPI_WTIME()

      sec=datestring2sec(DATEstring)
      call TimeInterpolation(sec,TC_CO2, BEFORE, AFTER, zweigh)

      iswap  = 0

! ----------------------- INITIALIZATION -------------
      IF (datestring.eq.DATESTART) then
!         scrivi qualche log sul numout

          CALL LOAD_CO2(TC_CO2%TimeStrings(TC_CO2%Before)) ! CALL trcrea(iperm1,bc)
          iswap = 1
          call swap_CO2


        CALL LOAD_CO2(TC_CO2%TimeStrings(TC_CO2%After))    ! CALL trcrea(iper,bc)


      ENDIF

! --------------------------------------------------------
! and now what we have to DO at every time step
! --------------------------------------------------------

! check the validity of the period in memory
      IF (BEFORE.ne.TC_CO2%Before) then
         TC_CO2%Before = BEFORE
         TC_CO2%After  = AFTER

         call swap_CO2
         iswap = 1


          CALL LOAD_CO2(TC_CO2%TimeStrings(TC_CO2%After))

          IF(lwp) WRITE (numout,*) 'Carbon dioxide factor DATA READ for Time = ', TC_CO2%TimeStrings(TC_CO2%After)
!      ******* LOADED NEW FRAME *************
      END IF

! compute the DATA at the given time step

      SELECT CASE (nsptint)
           CASE (0)  !  ------- no time interpolation
!      we have to initialize DATA IF we have changed the period
              IF (iswap.eq.1) THEN
                 zweigh = 1.0
                 call ACTUALIZE_CO2(zweigh)! initialize now fields with the NEW DATA READ
              END IF

          CASE (1) ! ------------linear interpolation ---------------
             call ACTUALIZE_CO2(zweigh)

      END SELECT


      bc_co2_partTime = MPI_WTIME()    - bc_co2_partTime
      bc_co2_TotTime  = bc_co2_TotTime + bc_co2_partTime

      END SUBROUTINE BC_CO2




! ******************************************************
!     SUBROUTINE LOAD_CO2(datestring)
!     loads BC/CO2_yyyy0107-00:00:00.nc in BC_mem.co2_dtatrc(:,2,:)
! ******************************************************
      SUBROUTINE LOAD_CO2(datestring)
         
         

          IMPLICIT NONE

          CHARACTER(LEN=17), INTENT(IN) :: datestring

!         local
          character(LEN=27) nomefile

          nomefile='BC/CO2_yyyy0107-00:00:00.nc'


    !     Starting I/O
    !    **********************************************************
         nomefile(1:27)='BC/CO2_'//datestring//'.nc'
         if(lwp) write(*,'(A,I4,A,A)') "LOAD_CO2 --> I am ", myrank, " starting reading forcing fields from ", nomefile(1:27)

         call readnc_slice_float_2d(nomefile,'CO2',buf2)
         co2_IO(:,:,2) = buf2*tmask(1,:,:)
      


      END SUBROUTINE LOAD_CO2


! ******************************************************
!     SUBROUTINE ACTUALIZE_CO2(zweigh)
!     performs time interpolation
!     x(1)*(1-zweigh) + x(2)*zweigh
! ******************************************************
      SUBROUTINE ACTUALIZE_CO2(zweigh)
        
         IMPLICIT NONE

         double precision, INTENT(IN) :: zweigh

        ! local
        INTEGER jj,ji

        DO ji=1,jpi
          DO jj=1,jpj
             ogstm_co2(jj,ji) = ( (1. - zweigh) * co2_IO(jj,ji,1)+ zweigh     * co2_IO(jj,ji,2) )
          END DO
        END DO


      END SUBROUTINE ACTUALIZE_CO2


! ****************************************************


      SUBROUTINE swap_CO2
         
         IMPLICIT NONE


         INTEGER jj,ji

            DO ji=1,jpi
              DO jj=1,jpj
                co2_IO(jj,ji,1) = co2_IO(jj,ji,2)
              END DO
            END DO


      END SUBROUTINE swap_CO2

END MODULE mod_cbc
