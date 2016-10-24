MODULE mod_gibbc

USE myalloc
USE BC_mem
USE TIME_MANAGER
USE calendar
USE mpi

implicit NONE

contains      
      
      SUBROUTINE BC_GIB(datestring)

      character(LEN=17), INTENT(IN) ::  datestring

! local declarations
! ==================
      REAL(8) sec,zweigh
      integer Before, After
      INTEGER iswap

      bc_gib_partTime = MPI_WTIME()

      sec=datestring2sec(DATEstring)
      call TimeInterpolation(sec,TC_GIB, BEFORE, AFTER, zweigh)

      iswap  = 0

! ----------------------- INITIALIZATION -------------
      IF (datestring.eq.DATESTART) then
!         scrivi qualche log sul numout

          CALL LOAD_GIB(TC_GIB%TimeStrings(TC_GIB%Before)) ! CALL trcrea(iperm1,bc)
          iswap = 1
          call swap_GIB


        CALL LOAD_GIB(TC_GIB%TimeStrings(TC_GIB%After))    ! CALL trcrea(iper,bc)


      ENDIF

! --------------------------------------------------------
! and now what we have to DO at every time step
! --------------------------------------------------------

! check the validity of the period in memory
      IF (BEFORE.ne.TC_GIB%Before) then
         TC_GIB%Before = BEFORE
         TC_GIB%After  = AFTER

         call swap_GIB
         iswap = 1


          CALL LOAD_GIB(TC_GIB%TimeStrings(TC_GIB%After))

          IF(lwp) WRITE (numout,*) 'Gibraltar factor DATA READ for Time = ', TC_GIB%TimeStrings(TC_GIB%After)

!      ******* LOADED NEW FRAME *************
      END IF

! compute the DATA at the given time step

      SELECT CASE (nsptint)
           CASE (0)  !  ------- no temporal interpolation
!      we have to initialize DATA IF we have changed the period
              IF (iswap.eq.1) THEN
                 zweigh = 1.0
                 call actualize_GIB(zweigh)! initialize now fields with the NEW DATA READ
              END IF

          CASE (1) ! ------------linear interpolation ---------------
             call actualize_GIB(zweigh)

      END SELECT


       bc_gib_partTime = MPI_WTIME()    - bc_gib_partTime
       bc_gib_TotTime  = bc_gib_TotTime + bc_gib_partTime


      END SUBROUTINE BC_GIB




! ******************************************************
!     SUBROUTINE LOAD_GIB(datestring)
!     loads BC/GIB_yyyy0107-00:00:00.nc in BC_mem.gib_dtatrc(:,2,:)
! ******************************************************
      SUBROUTINE LOAD_GIB(datestring)
         

          CHARACTER(LEN=17), INTENT(IN) :: datestring

!         local
          character(LEN=7)  nomevar
          character(LEN=27) nomefile
          INTEGER(4) jn,jv

          nomevar= '1234567'
          nomefile='BC/GIB_yyyy0107-00:00:00.nc'


    !     Starting I/O
    !    **********************************************************
          nomefile = 'BC/GIB_'//datestring//'.nc'
          if(lwp) write(*,'(A,I4,A,A)') "LOAD_GIB --> I am ", myrank, " starting reading forcing fields from ",nomefile(1:27)

          DO jn = 1, jn_gib
              nomevar = 'gib_'//ctrcnm(tra_matrix_gib(jn))
              if(lwp) write(*,*) "LOAD_GIB --> I am ", myrank,'name var is ', nomevar(1:7)
              !CALL ioogsnc_bc_1d2(nomefile,nomevar, Gsizeglo,gib_aux)
              CALL readnc_double_1d(nomefile,nomevar, Gsizeglo,gib_aux)

              DO jv=1,Gsize
                 gib_dtatrc(jv,2,jn) = gib_aux(gib_ridxt(1,jv))
              ENDDO
          ENDDO


      END SUBROUTINE LOAD_GIB

! ****************************************************

      SUBROUTINE actualize_GIB(zweigh)

         REAL(8), INTENT(IN) :: zweigh
!         local
         INTEGER jn, jv

         DO jn=1, jn_gib
             DO jv=1, Gsize
                 gib(jv,jn) = (1. - zweigh) * gib_dtatrc(jv,1,jn) + zweigh * gib_dtatrc(jv,2,jn)
             ENDDO
         ENDDO


      END SUBROUTINE actualize_GIB


! ****************************************************
      SUBROUTINE swap_GIB

!         local
          INTEGER jn, jv

          DO jn=1, jn_gib
              DO jv=1, Gsize
                  gib_dtatrc(jv,1,jn)=gib_dtatrc(jv,2,jn)
              ENDDO
          ENDDO

      END SUBROUTINE swap_GIB

END MODULE mod_gibbc