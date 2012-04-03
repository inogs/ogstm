! SUBROUTINE forcings_KEXT(datestring)
! si può copiare un header da dtadyn

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
      INTEGER iper,iperm1,iswap
      INTEGER jtask,kt



       dtadynparttime = MPI_WTIME()

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


! lb in any case, we need rhopn

!      jtask = 1; kt = 1 ; CALL eos(jtask,kt)


       dtadynparttime = MPI_WTIME() - dtadynparttime
       dtadyntottime = dtadyntottime + dtadynparttime



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

      IMPLICIT NONE

      CHARACTER(LEN=17), INTENT(IN) :: datestring

!      character*9  nomevar
      character*35 nomefile
!      INTEGER jinew, jjnew, jista, jiend, jjsta, jjend
!      INTEGER ji, jj, jk, ik

!      nomevar= '123456789'
      nomefile='FORCINGS/KextF_yyyy0107-00:00:00.nc'

!      nomevar(9:9)    = CHAR(0) ! once and for all
!      nomefile(36:36) = CHAR(0)


!     Starting I/O
!    **********************************************************
      nomefile(1:35) = 'FORCINGS/KextF_'//datestring//'.nc'

      if(lwp) write(*,*) "LOAD_KEXT --> I am ", rank, " starting reading forcing fields from ", nomefile(1:35)

      call readnc_slice_float(nomefile,'kextfact',buf2); kextIO(:,:,2) = buf2*tmask(:,:,1);
!      nomevar(1:8)='kextfact'; CALL ioogsnc_2dt2(nomefile,nomevar,jpiglo,jpjglo, zkef_f)
!     In this way only one file has to be read, also with more than a cpu
!      jista = nimpp
!      jiend = MIN(jpi+nimpp-1, jpiglo)
!      jjsta = njmpp
!      jjend = MIN(jpj+njmpp-1, jpjglo)
!
!
!
!      DO jj=jjsta, jjend
!         DO ji=jista, jiend
!            jinew = ji - nimpp + 1
!            jjnew = jj - njmpp + 1
!            kextIO(jinew,jjnew,2)=REAL(zkef_f(ji,jj),8)*tmask(jinew,jjnew,1)
!          END DO
!        END DO


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
        INTEGER ji,jj

        DO ji=1,jpi
          DO jj=1,jpj
             kef(ji,jj) = ( (1. - zweigh) * kextIO(ji,jj,1)+ zweigh     * kextIO(ji,jj,2) )
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


         INTEGER ji,jj

            DO ji=1,jpi
              DO jj=1,jpj
                kextIO(ji,jj,1) = kextIO(ji,jj,2)
              END DO
            END DO


      END SUBROUTINE swap_KEXT


