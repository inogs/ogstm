
      SUBROUTINE trcbio
!!!---------------------------------------------------------------------
!!!
!!!                       ROUTINE trcbio
!!!                     *******************
!!!
!!!  PURPOSE :
!!!  ---------
!!!     compute the now trend due to biogeochemical processes
!!!     and add it to the general trend of passive tracers equations.
!!!
!!!    Three options:
!!!
!!!   METHOD :
!!!   -------
!!!      each now biological flux is calculated  in FUNCTION of now
!!!      concentrations of tracers.
!!!      depending on the tracer, these fluxes are sources or sinks.
!!!      the total of the sources and sinks for each tracer
!!!      is added to the general trend.
!!!
!!!        tra = tra + zf...tra - zftra...
!!!                             |         |
!!!                             |         |
!!!                          source      sink
!!!
!!!
!!!      IF 'key_trc_diabio' key is activated, the biogeochemical
!!!    trends for passive tracers are saved for futher diagnostics.
!!!
!!!      multitasked on vertical slab (jj-loop)
!!!
!!!   MODIFICATIONS:
!!!   --------------

      USE myalloc
      USE myalloc_mpp
      USE BIO_mem
      USE BC_mem

       IMPLICIT NONE


!!!----------------------------------------------------------------------
!!! local declarations
!!! ==================
      LOGICAL sur,bot
      REAL(8) a(jptra),b(jptra),c(4),d(jptra_dia),er(10),d2(jptra_dia_2d)

      INTEGER jk,jj,ji,jb,jn
      INTEGER jtr,jtrmax,tra_idx

! omp variables
            INTEGER :: mytid, ntids

#ifdef __OPENMP1
            INTEGER ::  omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
            EXTERNAL :: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
#endif

!!!----------------------------------------------------------------------
!!! statement functions
!!! ===================


!   | --------------|
!   | BFM MODEL CALL|
!   | --------------|

        BIOparttime = MPI_WTIME()

#ifdef __OPENMP1
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
!!$omp    parallel default(none) private(mytid, jk,jj,ji) shared(sediPI,jpk,jpj,jpi,jn)
#ifdef __OPENMP1
        mytid = omp_get_thread_num()  ! take the thread ID
#endif
        IF (mytid+jn <= 4) then
            do jk=1,jpk
            do jj=1,jpj
            do ji=1,jpi
                sediPI(jk,jj,ji,jn+mytid)=0.
            end do
            end do
            end do
        ENDIF

!!$omp end parallel
         ENDDO

         DO jn=1, jptra_dia, ntids
!!$omp    parallel default(none) private(mytid, jk,jj,ji) shared(tra_DIA,jpk,jpj,jpi,jn)
#ifdef __OPENMP1
        mytid = omp_get_thread_num()  ! take the thread ID
#endif
        IF (mytid+jn <= jptra_dia) then
            do jk=1,jpk
            do jj=1,jpj
            do ji=1,jpi
               tra_DIA(jk,jj,ji,jn+mytid)=0.
            end do
            end do
            end do
        ENDIF
!!$omp end parallel
         ENDDO

      sediPI     = 0.
      tra_DIA    = 0.
      tra_DIA_2d = 0.



! $omp   parallel do default(none)  private(jb,jk,jj,ji,mytid,sur,bot,jtr,a,b,c,d,d2,er)
! $omp&      shared(NBFMPOINTS, BFMpoints,tra_idx,tra_matrix_gib,
! $omp&               restotr,jtrmax,trn,tn,sn,xpar,e3t,vatm,surf_mask,DAY_LENGTH,
! $omp&             sediPI,PH,tra_DIA,tra_DIA_2d,tra,rho,ice,co2,idxt2glo)

      MAIN_LOOP: DO  jb = 1, NBFMPOINTS
#ifdef __OPENMP1
        mytid = omp_get_thread_num()  ! take the thread ID
#endif

                 !IF( mytid + jb <= NBFMPOINTS ) THEN


                 ji = BFMpoints(1, jb)
                 jj = BFMpoints(2, jb)
                 jk = BFMpoints(3, jb)


                          sur = (jk .eq. 1)
                          bot = .FALSE.

                          DO jtr=1, jtrmax
                             a(jtr) = trn(jk,jj,ji,jtr) ! current biogeochemical concentrations
                          END DO
! Environmental regulating factors (er)

                          er(1)  = tn (jk,jj,ji)        ! Temperature (Celsius)
                          er(2)  = sn (jk,jj,ji)        ! Salinity PSU
                          er(3)  = rho(jk,jj,ji)        ! Density Kg/m3
                          er(4)  = ice                  ! from 0 to 1 adimensional
                          er(5)  = co2(jj,ji)           ! CO2 Mixing Ratios (ppm)  390
                          er(6)  = xpar(jk,jj,ji)       ! PAR umoles/m2/s | Watt to umoles photons W2E=1./0.217
                          er(7)  = DAY_LENGTH(jj,ji)    ! fotoperiod expressed in hours
                          er(8)  = e3t(jk,jj,ji)        ! depth in meters of the given cell
                          er(9)  = vatm(jj,ji) * surf_mask(jk) ! wind speed (m/s)
                          er(10) = PH(jk,jj,ji)         ! PH

                          call BFM0D_Input_EcologyDynamics(sur,bot,a,jtrmax,er)

                          call BFM0D_reset()

                         call EcologyDynamics()

                          if (sur) then
                             call BFM0D_Output_EcologyDynamics_surf(b, c, d ,d2)
                           else
                              call BFM0D_Output_EcologyDynamics(b, c, d)
                           endif

                          DO jtr=1, jtrmax
                             tra(jk,jj,ji,jtr) =tra(jk,jj,ji,jtr) +b(jtr) ! trend
                          END DO

                          DO jtr=1,4
                             sediPI(jk,jj,ji,jtr) = c(jtr) ! BFM output of sedimentation speed (m/d)
                          END DO

                          DO jtr=1,jptra_dia
                             tra_DIA(jk,jj,ji,jtr) = d(jtr) ! diagnostic
                          END DO

                          if (sur) then
                              DO jtr=1,jptra_dia_2d
                                 tra_DIA_2d(jj,ji,jtr) = d2(jtr) ! diagnostic
                              END DO
                          endif

                          PH(jk,jj,ji)=d(pppH) ! Follows solver guess, put 8.0 if pppH is not defined


             !ENDIF

                END DO MAIN_LOOP

! $omp end parallel do

                BIOparttime =  MPI_WTIME() -BIOparttime
                BIOtottime  = BIOtottime  + BIOparttime

      END SUBROUTINE trcbio
