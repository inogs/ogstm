
      SUBROUTINE trcbio
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE trcbio
CCC                     *******************
CCC
CCC  PURPOSE :
CCC  ---------
CCC     compute the now trend due to biogeochemical processes
CCC     and add it to the general trend of passive tracers equations.
CCC
CCC    Three options:
CCC
CC   METHOD :
CC   -------
CC      each now biological flux is calculated  in FUNCTION of now
CC      concentrations of tracers.
CC      depending on the tracer, these fluxes are sources or sinks.
CC      the total of the sources and sinks for each tracer
CC      is added to the general trend.
CC
CC        tra = tra + zf...tra - zftra...
CC                             |         |
CC                             |         |
CC                          source      sink
CC
CC
CC      IF 'key_trc_diabio' key is activated, the biogeochemical
CC    trends for passive tracers are saved for futher diagnostics.
CC
CC      multitasked on vertical slab (jj-loop)
CC
CCC   MODIFICATIONS:
CC   --------------

      USE myalloc
      USE myalloc_mpp
      USE BIO_mem
      USE BC_mem

       IMPLICIT NONE


CC----------------------------------------------------------------------
CC local declarations
CC ==================
      LOGICAL sur,bot
      REAL(8) a(jptra),b(jptra),c(4),d(jptra_dia-2),er(10),d2(jptra_dia_2d)

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
!$omp    parallel default(none) private(mytid, ji,jj,jk) shared(tra_DIA,jpi,jpj,jpk,jn)
#ifdef __OPENMP
        mytid = omp_get_thread_num()  ! take the thread ID
#endif
        IF (mytid+jn <= jptra_dia) then
            do jk=1,jpk
            do jj=1,jpj
            do ji=1,jpi
               tra_DIA(ji,jj,jk,jn+mytid)=0.
            end do
            end do
            end do
        ENDIF
!$omp end parallel
         ENDDO

      tra_DIA_2d=0.

      MAIN_LOOP: DO  jb = 1, NBFMPOINTS, ntids


!$omp   parallel default(none) private(ji,jj,jk,mytid,sur,bot,jtr,a,b,c,d,d2,er)
!$omp&      shared(jb,NBFMPOINTS, BFMpoints,tra_idx,tra_matrix_gib,
!$omp&               restotr,jtrmax,trn,tn,sn,xpar,e3t,vatm,surf_mask,DAY_LENGTH,
!$omp&             sediPI,PH,tra_DIA,tra_DIA_2d,tra,rho,ice,co2,idxt2glo,NPPF2)

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
                             a(jtr) = trn(ji,jj,jk,jtr) ! current biogeochemical concentrations
                          END DO
! Environmental regulating factors (er)

                          er(1)  = tn (ji,jj,jk)        ! Temperature (Celsius)
                          er(2)  = sn (ji,jj,jk)        ! Salinity PSU
                          er(3)  = rho(ji,jj,jk)        ! Density Kg/m3
                          er(4)  = ice                  ! from 0 to 1 adimensional
                          er(5)  = co2(ji,jj)           ! CO2 Mixing Ratios (ppm)  390
                          er(6)  = xpar(ji,jj,jk)       ! PAR umoles/m2/s | Watt to umoles photons W2E=1./0.217
                          er(7)  = DAY_LENGTH(ji,jj)    ! fotoperiod expressed in hours
                          er(8)  = e3t(jk)              ! depth in meters of the given cell
                          er(9)  = vatm(ji,jj) * surf_mask(jk) ! wind speed (m/s)
                          er(10) = PH(ji,jj,jk)         ! PH

                          call BFM0D_Input_EcologyDynamics(sur,bot,a,jtrmax,er)

                          call BFM0D_reset()

                          call EcologyDynamics()

                          if (sur) then
                             call BFM0D_Output_EcologyDynamics_surf(b, c, d ,d2)
                           else
                              call BFM0D_Output_EcologyDynamics(b, c, d)
                           endif

                          DO jtr=1, jtrmax
                             tra(ji,jj,jk,jtr) =tra(ji,jj,jk,jtr) +b(jtr) ! trend
                          END DO

                          DO jtr=1,4
                             sediPI(ji,jj,jk,jtr) = c(jtr) ! BFM output of sedimentation speed (m/d)
                          END DO

                          DO jtr=1,jptra_dia -2
                             tra_DIA(ji,jj,jk,jtr) = d(jtr) ! diagnostic
                          END DO

                          if (sur) then
                              DO jtr=1,jptra_dia_2d
                                 tra_DIA_2d(ji,jj,jtr) = d2(jtr) ! diagnostic
                              END DO
                          endif

                          PH(ji,jj,jk)=d(pppH) ! Follows solver guess, put 8.0 if pppH is not defined
                          NPPF2(ji,jj,jk)=d(ppF04) ! Flagellate production 



             ENDIF

!$omp end parallel

                END DO MAIN_LOOP

                BIOparttime =  MPI_WTIME() -BIOparttime
                BIOtottime  = BIOtottime  + BIOparttime

      END SUBROUTINE trcbio
