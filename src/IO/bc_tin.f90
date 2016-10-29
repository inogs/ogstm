MODULE mod_tinbc

USE myalloc
USE BC_mem
USE TIME_MANAGER
USE calendar
USE mpi

implicit NONE

contains      

      SUBROUTINE BC_TIN(datestring)

      character(LEN=17), INTENT(IN) ::  datestring

! local declarations
! ==================
      REAL(8) sec,zweigh
      integer Before, After
      INTEGER iswap

      bc_tin_partTime = MPI_WTIME()

      sec=datestring2sec(DATEstring)
      call TimeInterpolation(sec,TC_TIN, BEFORE, AFTER, zweigh)

      iswap  = 0

! ----------------------- INITIALIZATION -------------
      IF (datestring.eq.DATESTART) then

          CALL LOAD_TIN(TC_TIN%TimeStrings(TC_TIN%Before)) ! CALL trcrea(iperm1,bc)
          iswap = 1
          call swap_TIN


        CALL LOAD_TIN(TC_TIN%TimeStrings(TC_TIN%After))    ! CALL trcrea(iper,bc)


      ENDIF

! --------------------------------------------------------
! and now what we have to DO at every time step
! --------------------------------------------------------

! check the validity of the period in memory
      IF (BEFORE.ne.TC_TIN%Before) then
         TC_TIN%Before = BEFORE
         TC_TIN%After  = AFTER

         call swap_TIN
         iswap = 1


          CALL LOAD_TIN(TC_TIN%TimeStrings(TC_TIN%After))

          IF (lwp) WRITE (numout,*) 'Riverine factor DATA READ for Time = ', TC_TIN%TimeStrings(TC_TIN%After)

!      ******* LOADED NEW FRAME *************
      END IF

! compute the DATA at the given time step

      SELECT CASE (nsptint)
           CASE (0)  !  ------- no time interpolation
!      we have to initialize DATA IF we have changed the period
              IF (iswap.eq.1) THEN
                 zweigh = 1.0
                 call actualize_TIN(zweigh)! initialize now fields with the NEW DATA READ
              END IF

          CASE (1) ! ------------linear interpolation ---------------
             call actualize_TIN(zweigh)
      END SELECT



      bc_tin_partTime = MPI_WTIME()    - bc_tin_partTime
      bc_tin_TotTime  = bc_tin_TotTime + bc_tin_partTime


      END SUBROUTINE BC_TIN




! ******************************************************
!     SUBROUTINE LOAD_TIN(datestring)
!     loads BC/TIN_yyyy0107-00:00:00.nc in BC_mem.riv_dtatrc(:,2,:)
! ******************************************************
      SUBROUTINE LOAD_TIN(datestring)


          IMPLICIT NONE

          CHARACTER(LEN=17), INTENT(IN) :: datestring

!         local
          character(LEN=7)  nomevar
          character(LEN=27) nomefile
          INTEGER(4) jn,jv

          nomevar= '1234567'
          nomefile='BC/TIN_yyyy0107-00:00:00.nc'

    !     Starting I/O
    !    **********************************************************
         nomefile  ='BC/TIN_'//datestring//'.nc'
         if(lwp) write(*,'(A,I4,A,A)') "LOAD_TIN --> I am ", myrank, " starting reading forcing fields from ", nomefile(1:27)



          DO jn = 1, jn_riv

              nomevar = 'riv_'//ctrcnm(tra_matrix_riv(jn))
              !CALL ioogsnc_bc_1d2(nomefile,nomevar,Rsizeglo,riv_aux)
              CALL readnc_double_1d(nomefile,nomevar, Rsizeglo,riv_aux)
              DO jv=1,Rsize
                 riv_dtatrc(jv,2,jn) = riv_aux(riv_ridxt(1,jv))
              ENDDO
          ENDDO
          

      END SUBROUTINE LOAD_TIN

! ****************************************************

      SUBROUTINE actualize_TIN(zweigh)


         REAL(8), INTENT(IN) :: zweigh
!         local
         INTEGER jn, jv

         DO jn=1, jn_riv
             DO jv=1, Rsize
                 riv(jv,jn) = (1. - zweigh) * riv_dtatrc(jv,1,jn) + zweigh * riv_dtatrc(jv,2,jn)
             ENDDO
         ENDDO


      END SUBROUTINE actualize_TIN


! ****************************************************
      SUBROUTINE swap_TIN

!         local
          INTEGER jn, jv

          DO jn=1, jn_riv
              DO jv=1, Rsize
                  riv_dtatrc(jv,1,jn)=riv_dtatrc(jv,2,jn)
              ENDDO
          ENDDO

      END SUBROUTINE swap_TIN

END MODULE mod_tinbc
