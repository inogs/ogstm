      SUBROUTINE BC_GIB(datestring)


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


       dtatrcparttime = MPI_WTIME() - dtatrcparttime
       dtatrctottime = dtatrctottime + dtatrcparttime


      END SUBROUTINE BC_GIB




! ******************************************************
!     SUBROUTINE LOAD_GIB(datestring)
!     loads BC/GIB_yyyy0107-00:00:00.nc in BC_mem.gib_dtatrc(:,2,:)
! ******************************************************
      SUBROUTINE LOAD_GIB(datestring)
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
          INTEGER(4) jn,jv

          nomevar= '12345678'
          nomefile='BC/GIB_yyyy0107-00:00:00.nc'

          nomevar(8:8)    = CHAR(0) ! once and for all
          nomefile(28:28) = CHAR(0)

    !     Starting I/O
    !    **********************************************************
          nomefile(1:27)='BC/GIB_'//datestring//'.nc'
          if(lwp) write(*,*) "LOAD_GIB --> I am ", rank, " starting reading forcing fields from ",nomefile(1:27)





          DO jn = 1, jn_gib

              nomevar(1:7) = 'gib_'//ctrcnm(tra_matrix_gib(jn))(1:3)
              if(lwp) write(*,*) "LOAD_GIB --> I am ", rank,'name var is ', nomevar(1:7)
              CALL ioogsnc_bc_1d2(nomefile,nomevar, Gsizeglo,gib_aux)

              DO jv=1,Gsize
                 gib_dtatrc(jv,2,jn) = gib_aux(gib_ridxt(jv,1))
              ENDDO
          ENDDO


      END SUBROUTINE LOAD_GIB

! ****************************************************

      SUBROUTINE actualize_GIB(zweigh)
         USE myalloc
         USE BC_mem
         IMPLICIT NONE

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
          use myalloc ! oink oink
          use BC_mem  ! for jn_gib, Gsize, gib_dtatrc
          IMPLICIT NONE

!         local
          INTEGER jn, jv

          DO jn=1, jn_gib
              DO jv=1, Gsize
                  gib_dtatrc(jv,1,jn)=gib_dtatrc(jv,2,jn)
              ENDDO
          ENDDO

      END SUBROUTINE swap_GIB
