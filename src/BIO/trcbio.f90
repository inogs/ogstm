
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
      ! epascolo USE myalloc_mpp
      USE BIO_mem
      USE BC_mem
      USE mpi
      
      IMPLICIT NONE


!!!----------------------------------------------------------------------
!!! local declarations
!!! ==================
      logical :: sur,bot
      double precision,dimension(jptra) :: a,b
      double precision,dimension(4) :: c
      double precision,dimension(jptra_dia) :: d
      double precision,dimension(10) :: er
      double precision,dimension(jptra_dia_2d) :: d2

      integer :: jk,jj,ji,jb,jn
      integer :: jtr,jtrmax,tra_idx


!!!----------------------------------------------------------------------
!!! statement functions
!!! ===================


!   | --------------|
!   | BFM MODEL CALL|
!   | --------------|

        BIOparttime = MPI_WTIME()

          surf_mask(:) = 0.
          surf_mask(1) = 1.
! -------------------------------------------------

          tra_idx = tra_matrix_gib(1)
          jtrmax=jptra

! ---------------- Fuori dai punti BFM

      sediPI     = 0.
      tra_DIA    = 0.
      tra_DIA_2d = 0.

! $omp   parallel do default(none)  private(jb,jk,jj,ji,mytid,sur,bot,jtr,a,b,c,d,d2,er)
! $omp&      shared(NBFMPOINTS, BFMpoints,tra_idx,tra_matrix_gib,
! $omp&               restotr,jtrmax,trn,tn,sn,xpar,e3t,vatm,surf_mask,DAY_LENGTH,
! $omp&             sediPI,PH,tra_DIA,tra_DIA_2d,tra,rho,ice,co2,idxt2glo)

      MAIN_LOOP: DO  jb = 1, NBFMPOINTS

                 !IF( mytid + jb <= NBFMPOINTS ) THEN


                 ji = BFMpoints(3, jb)
                 jj = BFMpoints(2, jb)
                 jk = BFMpoints(1, jb)


                          sur = (jk .eq. 1)
                          bot = .FALSE.

                          DO jtr=1, jtrmax
                             a(jtr) = trn(jk,jj,ji,jtr) ! current biogeochemical concentrations
                        !      WRITE(*,200) ,'I',jk,jj,ji,jtr,trn(jk,jj,ji,jtr)
                             
                          END DO
! Environmental regulating factors (er)

                          er(1)  = tn (jk,jj,ji)        ! Temperature (Celsius)
                          er(2)  = sn (jk,jj,ji)        ! Salinity PSU
                          er(3)  = rho(jk,jj,ji)        ! Density Kg/m3
                          er(4)  = ice                  ! from 0 to 1 adimensional
                          er(5)  = ogstm_co2(jj,ji)           ! CO2 Mixing Ratios (ppm)  390
                          er(6)  = xpar(jk,jj,ji)       ! PAR umoles/m2/s | Watt to umoles photons W2E=1./0.217
                          er(7)  = DAY_LENGTH(jj,ji)    ! fotoperiod expressed in hours
                          er(8)  = e3t(jk,jj,ji)        ! depth in meters of the given cell
                          er(9)  = vatm(jj,ji) * surf_mask(jk) ! wind speed (m/s)
                          er(10) = ogstm_PH(jk,jj,ji)         ! PH
                        !   WRITE(*,201),'ERR',jk,jj,ji,er(:)
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
                        !      WRITE(*,200),'OB',jk,jj,ji,jtr,b(jtr)
                        !      WRITE(*,200),'TRA',jk,jj,ji,jtr,tra(jk,jj,ji,jtr)
                          END DO

                          DO jtr=1,4
                             ogstm_sediPI(jk,jj,ji,jtr) = c(jtr) ! BFM output of sedimentation speed (m/d)
                          END DO

                          DO jtr=1,jptra_dia
                             tra_DIA(jtr,jk,jj,ji) = d(jtr) ! diagnostic
                        !      WRITE(*,200),'IO3',jk,jj,ji,jtr,tra_DIA(jk,jj,ji,jtr)
                          END DO
                          if (sur) then
                              DO jtr=1,jptra_dia_2d
                                 tra_DIA_2d(jj,ji,jtr) = d2(jtr) ! diagnostic
                              !    WRITE(*,202),'IO2',jj,ji,jtr,tra_DIA_2d(jj,ji,jtr)
                              END DO
                          endif

                          ogstm_PH(jk,jj,ji)=d(pppH) ! Follows solver guess, put 8.0 if pppH is not defined


             !ENDIF

                END DO MAIN_LOOP
! 200    FORMAT(' ',A3,I4,I4,I4,I4,D30.23)
! 201    FORMAT(' ',A3,I4,I4,I4,10D30.23)
! 202    FORMAT(' ',A3,I4,I4,I4,10D30.23)
                
                          
! $omp end parallel do

                BIOparttime =  MPI_WTIME() -BIOparttime
                BIOtottime  = BIOtottime  + BIOparttime
               
      END SUBROUTINE trcbio
