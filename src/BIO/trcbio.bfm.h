
       USE myalloc
       USE myalloc_mpp
       USE BIO_mem
       USE BC_mem

       IMPLICIT NONE


CC----------------------------------------------------------------------
CC local declarations
CC ==================
      LOGICAL sur,bot
      REAL(8) a(jptra),b(jptra),c(4),d(jptra_dia),er(10)

      INTEGER ji,jj,jk,jb,jn
      INTEGER jtr,jtrmax,tra_idx

! omp variables
            INTEGER :: mytid, ntids

#ifdef __OPENMP
            INTEGER ::  omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
            EXTERNAL :: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
#endif

CC----------------------------------------------------------------------
CC statement functions
CC ===================


C   | --------------|
C   | BFM MODEL CALL|
C   | --------------|
C
        BIOparttime = MPI_WTIME()

#ifdef __OPENMP
      ntids = omp_get_max_threads() ! take the number of threads
      mytid = -1000000
#else
      ntids = 1
      mytid = 0
#endif
          surf_mask(:) = 0.
          surf_mask(1) = 1.
! -------------------------------------------------

          tra_idx = tra_matrix_gib(1)
          jtrmax=jptra

! ---------------- Fuori dai punti BFM
         DO jn=1,4,ntids
!$omp    parallel default(none) private(mytid, ji,jj,jk) shared(sediPI,jpi,jpj,jpk,jn)
#ifdef __OPENMP
        mytid = omp_get_thread_num()  ! take the thread ID
#endif
        IF (mytid+jn <= 4) then
            do jk=1,jpk
            do jj=1,jpj
            do ji=1,jpi
                sediPI(ji,jj,jk,jn+mytid)=0.
            end do
            end do
            end do
        ENDIF

!$omp end parallel
         ENDDO

         DO jn=1, jptra_dia, ntids
!$omp    parallel default(none) private(mytid, ji,jj,jk) shared(tra_pp,jpi,jpj,jpk,jn)
#ifdef __OPENMP
        mytid = omp_get_thread_num()  ! take the thread ID
#endif
        IF (mytid+jn <= jptra_dia) then
            do jk=1,jpk
            do jj=1,jpj
            do ji=1,jpi
               tra_pp(ji,jj,jk,jn+mytid)=0.
            end do
            end do
            end do
        ENDIF
!$omp end parallel
         ENDDO



!      MAIN_LOOP: DO  jj = 1, jpjm1, ntids
      MAIN_LOOP: DO  jb = 1, NBFMPOINTS, ntids


!$omp   parallel default(none) private(ji,jj,jk,mytid,sur,bot,jtr,a,b,c,d,er)
!$omp&      shared(jb,NBFMPOINTS, BFMpoints,tra_idx,tra_matrix_gib,
!$omp&               restotr,jtrmax,trn,tn,sn,xpar,e3t,vatm,surf_mask,DAY_LENGTH,
!$omp&             sediPI,PH,tra_pp,tra,rho,opa_ice,opa_co2,idxt2glo)

#ifdef __OPENMP
        mytid = omp_get_thread_num()  ! take the thread ID
#endif

                 IF( mytid + jb <= NBFMPOINTS ) THEN


                 ji = BFMpoints(1, mytid+jb)
                 jj = BFMpoints(2, mytid+jb)
                 jk = BFMpoints(3, mytid+jb)


                          sur = (jk .eq. 1)
                          bot = .FALSE.

                          DO jtr=1, jtrmax
                             a(jtr) = trn(ji,jj,jk,jtr)
                          END DO
! Environmental regulating factors (er)
                          er(1)  = tn (ji,jj,jk)
                          er(2)  = sn (ji,jj,jk)
                          er(3)  = rho(ji,jj,jk)
                          er(4)  = opa_ice
                          er(5)  = opa_co2(ji,jj)
                          er(6)  = xpar(ji,jj,jk)
                          er(7)  = DAY_LENGTH(ji,jj)
                          er(8)  = e3t(jk)
                          er(9)  = vatm(ji,jj) * surf_mask(jk)
                          er(10) = PH(ji,jj,jk)
                          call BFM0D_Input_EcologyDynamics(sur,bot,a,jtrmax,er)

                          call BFM0D_reset()

                          call EcologyDynamics()

                          call BFM0D_Output_EcologyDynamics(b, c, d)

                          DO jtr=1, jtrmax
                             tra(ji,jj,jk,jtr) =tra(ji,jj,jk,jtr) +b(jtr) ! trend
                          END DO

                          DO jtr=1,4
                             sediPI(ji,jj,jk,jtr) = c(jtr) ! sedimentation velocities
                          END DO

! Last record (jptra_dia) for evaporation rates

                          DO jtr=1,jptra_dia
                             tra_pp(ji,jj,jk,jtr) = d(jtr) ! diagnostic
                          END DO

                          PH(ji,jj,jk)=d(pppH) ! Follows solver guess, put 8.0 if pppH is not defined


             ENDIF

!$omp end parallel

                END DO MAIN_LOOP

                BIOparttime =  MPI_WTIME() -BIOparttime
                BIOtottime  = BIOtottime  + BIOparttime
