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

      dtadynparttime = MPI_WTIME()

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


       dtatrcparttime = MPI_WTIME() - dtatrcparttime
       dtatrctottime = dtatrctottime + dtatrcparttime


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
          character*8  nomevar
          character*28 nomefile
          INTEGER tra_matrix(51)
          INTEGER(4) jn,jv

          nomevar= '123456789'
          nomefile='BC/ATM_yyyy0107-00:00:00.nc'

          nomevar(8:8)    = CHAR(0) ! once and for all
          nomefile(28:28) = CHAR(0)

    !     Starting I/O
    !    **********************************************************
         nomefile(1:27)='BC/ATM_'//datestring//'.nc'
         if(lwp) write(*,*) "LOAD_ATM --> I am ", rank, " starting reading forcing fields from ", nomefile(1:27)




          tra_matrix(1) = ppO2o
          tra_matrix(2) = ppN1p
          tra_matrix(3) = ppN3n
          tra_matrix(4) = ppN5s

          DO jn = 1, jn_atm ! 4

              nomevar(1:7) = 'atm_'//ctrcnm(tra_matrix(jn))(1:3)
              CALL ioogsnc_bc_1d2(nomefile,nomevar,Asizeglo,atm_aux)
              DO jv=1,Asize
                 atm_dtatrc(jv,2,jn) = atm_aux(atm_ridxt(jv,1))
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
