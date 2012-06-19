
       USE myalloc
       USE myalloc_mpp
       USE TIME_MANAGER
       USE BIO_mem
       USE BC_mem

       IMPLICIT NONE


CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER kt,ktask
      LOGICAL isDumpAscii,sur,bot,isT,isBIO
      REAL(8) a(jptra),b(jptra),c(4),d(jptra_dia),er(10)

      INTEGER ji,jj,jk!,jn
      INTEGER gji,gjj,gjk
      INTEGER CSi,CS(10,3)
      INTEGER jtr,jtrmax,tra_idx
! omp variables
            INTEGER :: mytid, ntids!, itid

#ifdef __OPENMP
            INTEGER ::  omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
            EXTERNAL :: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
#endif

CC----------------------------------------------------------------------
CC statement functions
CC ===================


      INTERFACE OPA_Output_EcologyDynamics
         subroutine OPA_Output_EcologyDynamics(opa_tra, dim_opa_tra, sediPI, local_opa_dia)
!            use global_mem, ONLY:RLEN
            IMPLICIT NONE
            integer dim_opa_tra
            real(8):: sediPI(4)
            real(8):: opa_tra(dim_opa_tra)
            real(8):: local_opa_dia(33)
         end subroutine
      END INTERFACE

C   | --------------|
C   | BFM MODEL CALL|
C   | --------------|
C
#ifdef __OPENMP
      ntids = omp_get_max_threads() ! take the number of threads
      mytid = -1000000
#else
      ntids = 1
      mytid = 0
#endif
          surf_mask(:) = 0.
          surf_mask(1) = 1.


          tra_idx = tra_matrix_gib(1)
          jtrmax=jptra

          isDumpAscii=IsAnAveDump(NOW_datestring)

      MAIN_LOOP: DO  jj = 1, jpjm1, ntids

!$omp   parallel default(none) private(jk,ji,mytid,isT,isBIO,sur,bot,jtr,a,b,c,d,er,gji,gjj,gjk,CSi,CS)
!$omp&      shared(jj,jpjm1,jpkbm1,jpim1,Tmask,tra_idx,tra_matrix_gib,
!$omp&               restotr,jtrmax,trn,tn,sn,xpar,e3t,vatm,surf_mask,DAY_LENGTH,
!$omp&             sediPI,PH,tra_pp,tra,rhopn,opa_ice,opa_co2,idxt2glo,isDumpAscii,NOW_datestring)

#ifdef __OPENMP
        mytid = omp_get_thread_num()  ! take the thread ID
#endif

                 IF( mytid + jj <= jpjm1 ) THEN

                 DO jk=1,jpkbm1
                    DO ji = 2,jpim1

                       isT     = Tmask(ji,jj+mytid,jk) .eq. 1
                       isBIO   = restotr(ji,jj+mytid,jk,tra_idx) .eq. 0 ! no nudging points are considered

                       IF( isT .and. isBIO) THEN

                          sur = (jk .eq. 1)

                          bot = .FALSE.

                          DO jtr=1, jtrmax
                             a(jtr) = trn(ji,jj+mytid,jk,jtr)
                          END DO
! Environmental regulating factors (er)
                          er(1)  = tn(ji,jj+mytid,jk)
                          er(2)  = sn(ji,jj+mytid,jk)
                          er(3)  = rhopn(ji,jj+mytid,jk)
                          er(4)  = opa_ice
                          er(5)  = opa_co2
                          er(6)  = xpar(ji,jj+mytid,jk)
                          er(7)  = DAY_LENGTH(ji,jj+mytid)
                          er(8)  = e3t(jk)
                          er(9)  = vatm(ji,jj+mytid) * surf_mask(jk)
                          er(10) = PH(ji,jj+mytid,jk)

                          call OPA_Input_EcologyDynamics(sur,bot,a,jtrmax,er)

                          call OPA_reset()

                          call EcologyDynamics()

                          call OPA_Output_EcologyDynamics(b, jtrmax, c, d)

                          DO jtr=1, jtrmax
                             tra(ji,jj+mytid,jk,jtr) =tra(ji,jj+mytid,jk,jtr) +b(jtr) ! trend
                          END DO

                          DO jtr=1,4
                             sediPI(ji,jj+mytid,jk,jtr) = c(jtr) ! sedimentation velocities
                          END DO

                          DO jtr=1,jptra_dia
                             tra_pp(ji,jj+mytid,jk,jtr) = d(jtr) ! diagnostic
                          END DO

                          PH(ji,jj+mytid,jk)=d(9) ! Follows solver guess

! Diagnostic
                          IF (isDumpAscii) THEN

! From local to global indexing

                             gji = idxt2glo(ji, jj, jk,1)
                             gjj = idxt2glo(ji, jj, jk,2)
                             gjk = idxt2glo(ji, jj, jk,3)

! CheckSites (CS) localization --- CS(8,3)

                             CS(1,1)=36;  CS(1,2)=46    ! Alboran Sea

                             CS(2,1)=120; CS(2,2)=80    ! West Med

                             CS(3,1)=170; CS(3,2)=80    ! Tyrrhenian

                             CS(4,1)=131; CS(4,2)=105   ! DYFAMED

                             CS(5,1)=220; CS(5,2)=45    ! Ionian

                             CS(6,1)=310; CS(6,2)=30    ! East Med

                             CS(7,1)=172; CS(7,2)=122   ! North Adriatic

                             CS(8,1)=212; CS(8,2)=93    ! South Adriatic

                             CS(9,1)=266; CS(9,2)=78    ! North Aegean

                             CS(10,1)=279; CS(10,2)=81  ! North Aegean Critical

                             DO CSi=1,10
                                IF ( (gji .EQ. CS(CSi,1)) .AND. (gjj .EQ. CS(CSi,2)) ) THEN
                                   CALL OPA_SS_OUTPUT(gji,gjj,gjk,NOW_datestring) ! BFM routine that dumps all key fluxes
                                ENDIF
                             ENDDO
! Plot BFM status when carbonate system fails to converge PH_ = -1
                             IF ( d(9) .LT. 0. ) THEN
                                CALL OPA_SS_OUTPUT(gji,gjj,-1,NOW_datestring) ! BFM routine that dumps all key fluxes
                             ENDIF
                          ENDIF
! End Diagnostic

                       ELSE
                          sediPI(ji,jj+mytid,jk,:)=0
                          tra_pp(ji,jj+mytid,jk,:)=0
                       ENDIF

                    END DO

                END DO

             ENDIF

!$omp end parallel

                END DO MAIN_LOOP

