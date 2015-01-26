      SUBROUTINE BC_ATM(datestring)


       USE myalloc
       USE myalloc_mpp
       USE BC_mem
       USE TIME_MANAGER
       IMPLICIT NONE

      character(LEN=17), INTENT(IN) ::  datestring

! local declarations
! ==================
      REAL(8) sec,zweigh
      integer Before, After
      INTEGER iswap

      bc_atm_partTime = MPI_WTIME()

      sec=datestring2sec(DATEstring)
      call TimeInterpolation(sec,TC_ATM, BEFORE, AFTER, zweigh)

      iswap  = 0

! ----------------------- INITIALIZATION -------------
      IF (datestring.eq.DATESTART) then
!         scrivi qualche log sul numout

          CALL LOAD_ATM(TC_ATM%TimeStrings(TC_ATM%Before)) ! CALL trcrea(iperm1,bc)
          iswap = 1
          call swap_ATM


        CALL LOAD_ATM(TC_ATM%TimeStrings(TC_ATM%After))    ! CALL trcrea(iper,bc)


      ENDIF

! --------------------------------------------------------
! and now what we have to DO at every time step
! --------------------------------------------------------

! check the validity of the period in memory
      IF (BEFORE.ne.TC_ATM%Before) then
         TC_ATM%Before = BEFORE
         TC_ATM%After  = AFTER

         call swap_ATM
         iswap = 1


          CALL LOAD_ATM(TC_ATM%TimeStrings(TC_ATM%After))

          IF(lwp) WRITE (numout,*) 'Atmospheric factor DATA READ for Time = ', TC_ATM%TimeStrings(TC_ATM%After)
!      ******* LOADED NEW FRAME *************
      END IF

! compute the DATA at the given time step

      SELECT CASE (nsptint)
           CASE (0)  !  ------- no time interpolation
!      we have to initialize DATA IF we have changed the period
              IF (iswap.eq.1) THEN
                 zweigh = 1.0
                 call ACTUALIZE_ATM(zweigh)! initialize now fields with the NEW DATA READ
              END IF

          CASE (1) ! ------------linear interpolation ---------------
             call ACTUALIZE_ATM(zweigh)

      END SELECT


      bc_atm_partTime = MPI_WTIME()    - bc_atm_partTime
      bc_atm_TotTime  = bc_atm_TotTime + bc_atm_partTime

      END SUBROUTINE BC_ATM




! ******************************************************
!     SUBROUTINE LOAD_ATM(datestring)
!     loads BC/ATM_yyyy0107-00:00:00.nc in BC_mem.atm_dtatrc(:,2,:)
! ******************************************************
      SUBROUTINE LOAD_ATM(datestring)
          USE calendar
          USE myalloc
          USE myalloc_mpp
          USE BC_mem
          USE TIME_MANAGER

          IMPLICIT NONE

          CHARACTER(LEN=17), INTENT(IN) :: datestring

!         local
          character(LEN=7)  nomevar
          character(LEN=27) nomefile
          INTEGER(4) jn,jv

          nomevar= '1234567'
          nomefile='BC/ATM_yyyy0107-00:00:00.nc'


    !     Starting I/O
    !    **********************************************************
         nomefile(1:27)='BC/ATM_'//datestring//'.nc'
         if(lwp) write(*,'(A,I4,A,A)') "LOAD_ATM --> I am ", rank, " starting reading forcing fields from ", nomefile(1:27)





          DO jn = 1, jn_atm

              nomevar = 'atm_'//ctrcnm(tra_matrix_atm(jn))
              !CALL ioogsnc_bc_1d2(nomefile,nomevar,Asizeglo,atm_aux)
              CALL readnc_double_1d(nomefile,nomevar, Asizeglo,atm_aux)
              DO jv=1,Asize
                 atm_dtatrc(jv,2,jn) = atm_aux(atm_ridxt(1,jv))
              ENDDO
          ENDDO


      END SUBROUTINE LOAD_ATM

! ****************************************************

      SUBROUTINE actualize_ATM(zweigh)
         USE myalloc
         USE BC_mem
         IMPLICIT NONE

         REAL(8), INTENT(IN) :: zweigh
!         local
         INTEGER jn, jv

         DO jn=1, jn_atm
             DO jv=1, Asize
                 atm(jv,jn) = (1. - zweigh) * atm_dtatrc(jv,1,jn) + zweigh * atm_dtatrc(jv,2,jn)
             ENDDO
         ENDDO


      END SUBROUTINE actualize_ATM


! ****************************************************
      SUBROUTINE swap_ATM
          use myalloc ! oink oink
          use BC_mem  ! for jn_atm, Asize, atm_dtatrc
          IMPLICIT NONE

!         local
          INTEGER jn, jv

          DO jn=1, jn_atm
              DO jv=1, Asize
                  atm_dtatrc(jv,1,jn)=atm_dtatrc(jv,2,jn)
              ENDDO
          ENDDO

      END SUBROUTINE swap_ATM
