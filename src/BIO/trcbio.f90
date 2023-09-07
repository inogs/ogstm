
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
      USE BIO_mem
      USE OPT_mem, ONLY: PAR,SWR_RT
      USE BC_mem
      USE mpi

! ----------------------------------------------------------------------
!  BEGIN BC_REFACTORING SECTION
!  ---------------------------------------------------------------------

      use bc_set_mod

! ----------------------------------------------------------------------
!  END BC_REFACTORING SECTION
!  ---------------------------------------------------------------------
      
      IMPLICIT NONE


!!!----------------------------------------------------------------------
!!! local declarations
!!! ==================

#ifdef BFMv2
      logical :: sur,bot
      double precision,dimension(jptra) :: a,b
      double precision,dimension(4) :: c
      double precision,dimension(jptra_dia) :: d
      double precision,dimension(10) :: er
      double precision,dimension(jptra_dia_2d) :: d2
#else
      double precision,dimension(jptra,jpk) :: b
      double precision,dimension(jpk,jptra) :: a
      double precision,dimension(4,jpk) :: c
      double precision,dimension(jptra_dia,jpk) :: d
      double precision,dimension(jpk,18) :: er
      double precision,dimension(jptra_dia_2d) :: d2
#endif


      integer :: jk,jj,ji,jb,jn
      integer :: jtr,jtrmax,tra_idx
      integer :: bottom
      double precision :: correct_fact


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

          ! tra_idx = tra_matrix_gib(1)
          jtrmax=jptra

! ---------------- Fuori dai punti BFM

      ogstm_sediPI=0.
      tra_DIA    = 0.
      tra_DIA_2d = 0. ! da sistemare

#ifdef BFMv2
      tra_DIA_2d = 0.



      MAIN_LOOP: DO  jb = 1, NBFMPOINTS


                 ji = BFMpoints(3, jb)
                 jj = BFMpoints(2, jb)
                 jk = BFMpoints(1, jb)


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
                          er(5)  = ogstm_co2(jj,ji)           ! CO2 Mixing Ratios (ppm)  390
                          if (is_night(COMMON_DATEstring)) then
                              er(6)  = 0.001       ! PAR umoles/m2/s | Watt to umoles photons W2E=1./0.217
                          else
                              er(6)  = 2.0 * xpar(jk,jj,ji)       ! PAR umoles/m2/s | Watt to umoles photons W2E=1./0.217
                          endif
                          er(7)  = DAY_LENGTH(jj,ji)    ! fotoperiod expressed in hours
                          er(8)  = e3t(jk,jj,ji)        ! depth in meters of the given cell
                          er(9)  = vatm(jj,ji) * surf_mask(jk)  ! wind speed (m/s)
                          er(10) = ogstm_PH(jk,jj,ji)         ! PH

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
                             ogstm_sediPI(jk,jj,ji,jtr) = c(jtr) ! BFM output of sedimentation speed (m/d)
                          END DO

                          DO jtr=1,jptra_dia -2 ! We skip the last two ppHT1 and ppHT2
                             tra_DIA(jtr,jk,jj,ji) = d(jtr) ! diagnostic
                          END DO

                          if (sur) tra_DIA_2d(:,jj,ji) = d2(:) ! diagnostic

                          ogstm_PH(jk,jj,ji)=d(pppH) ! Follows solver guess, put 8.0 if pppH is not defined

                          NPPF2(jk,jj,ji)=d(ppF04) ! Flagellate production

                END DO MAIN_LOOP
#else

!    Initialization
      a        = 1.0
      er       = 1.0
      er(:,10) = 8.1

#ifdef gdept1d
! er(:,11) is calculated outside the loop on ji,jj
      do jk=1, jpk
         correct_fact= 1.0D0

         if ( (gdept(jk) .GT. 1000.0D0 ) .AND.  (gdept(jk) .LT. 2000.0D0 )) then
             correct_fact= 0.25D0
         endif

         if (gdept(jk) .GE. 2000.0D0 ) then
             correct_fact= 0.0D0
         endif

         er(jk,16) = correct_fact * ( gdept(jpk)-gdept(jk) ) /gdept(jpk)
      enddo
#endif


      DO ji=1,jpi
      DO jj=1,jpj
      if (bfmmask(1,jj,ji) == 0) CYCLE
      bottom = mbathy(jj,ji)


                          DO jtr=1, jtrmax

                             a(1:bottom, jtr) = trn(1:bottom,jj,ji,jtr) ! current biogeochemical concentrations

                          END DO

! Environmental regulating factors (er,:)

                          er(1:bottom,1)  = tn (1:bottom,jj,ji)   ! Temperature (Celsius)
                          er(1:bottom,2)  = sn (1:bottom,jj,ji)   ! Salinity PSU
                          er(1:bottom,3)  = rho(1:bottom,jj,ji)   ! Density Kg/m3
                          er(1       ,4)  = ice                   ! from 0 to 1 adimensional
                          er(1       ,5)  = ogstm_co2(jj,ji)      ! CO2 Mixing Ratios (ppm)  390
                          if (is_night(COMMON_DATEstring)) then
!                             er(6)  = 0.001       ! PAR umoles/m2/s | Watt to umoles photons W2E=1./0.217
                              er(1:bottom,6)  = 0.001 ! PAR for diatoms
                              er(1:bottom,7)  = 0.001 ! PAR for flagellates
                              er(1:bottom,8)  = 0.001 ! PAR for pico phytoplankton
                              er(1:bottom,9)  = 0.001 ! PAR for dinoflagellates
                              er(1:bottom,10) = 0.001 ! total PAR for CDOM
                              er(1:bottom,11) = 0.001 * 0.217/0.5!  Short wave radiation [0 4 um] W/m2
                          else
!                             er(6)  = 2.0 * xpar(jk,jj,ji)       ! PAR umoles/m2/s | Watt to umoles photons W2E=1./0.217
                              do jk=1, bottom
                                  er(jk,6)  = instant_par(COMMON_DATEstring,xpar(jk,jj,ji))           ! PAR for diatoms
                                  er(jk,7)  = instant_par(COMMON_DATEstring,xpar(jk,jj,ji))           ! PAR for flagellates
                                  er(jk,8)  = instant_par(COMMON_DATEstring,xpar(jk,jj,ji))           ! PAR for pico phytoplankton
                                  er(jk,9)  = instant_par(COMMON_DATEstring,xpar(jk,jj,ji))           ! PAR for dinoflagellates
                                  er(jk,10) = instant_par(COMMON_DATEstring,xpar(jk,jj,ji))           ! total PAR for CDOM
                                  er(jk,11) = instant_par(COMMON_DATEstring,xpar(jk,jj,ji))*0.217/0.5 !  Short wave radiation [0 4 um] W/m2
                              enddo
                          endif
!                         er(jk,6) = instant_par(COMMON_DATEstring,xpar(jk,jj,ji))  ! PAR umoles/m2/s | Watt to umoles photons W2E=1./0.217
                          er(1       ,12)  = DAY_LENGTH(jj,ji)    ! fotoperiod expressed in hours
                          er(1:bottom,13)  = e3t(1:bottom,jj,ji)        ! depth in meters of the given cell
                          er(1       ,14)  = vatm(jj,ji)                ! wind speed (m/s)
                          er(1:bottom,15) = ogstm_PH(1:bottom,jj,ji)   ! 8.1
#ifndef gdept1d
                         do jk=1,bottom
                             correct_fact= 1.0D0
                             if ( (gdept(jk,jj,ji) .GT. 1000.0D0 ) .AND.  (gdept(jk,jj,ji) .LT. 2000.0D0)) then
                                 correct_fact= 0.25D0
                             endif

                             if (gdept(jk,jj,ji) .GE. 2000.0D0 ) then
                                 correct_fact= 0.0D0
                             endif

                             er(jk,16) = correct_fact * ( gdept(jpk,jj,ji)-gdept(jk,jj,ji) ) /gdept(jpk,jj,ji)
                         enddo
#endif
                         er(1       ,17)  = atm_Hg0(jj,ji)                   ! ATM Hg0    
                         er(1       ,18)  = atm_WIN(jj,ji)               !wind from atm model

                          call BFM1D_Input_EcologyDynamics(bottom,a,jtrmax,er)

                         call BFM1D_reset()

                         call EcologyDynamics()

                         call BFM1D_Output_EcologyDynamics(b, c, d, d2)

                          DO jtr=1, jtrmax
                             tra(1:bottom,jj,ji,jtr) =tra(1:bottom,jj,ji,jtr) +b(jtr,1:bottom) ! trend
                          END DO

                          DO jtr=1,4
                             ogstm_sediPI(1:bottom,jj,ji,jtr) = c(jtr,1:bottom)      ! BFM output of sedimentation speed (m/d)
                          END DO


                          DO jk = 1,bottom
                          DO jtr=1,jptra_dia
                             tra_DIA(jtr, jk ,jj,ji) = d(jtr,jk) ! diagnostic
                          END DO
                          ENDDO

                         tra_DIA_2d(:,jj,ji) = d2(:) ! diagnostic

                          ogstm_PH(1:bottom,jj,ji) = d(pppH,1:bottom) ! Follows solver guess, put 8.0 if pppH is not defined

      END DO
      END DO

! ----------------------------------------------------------------------
!  BEGIN BC_REFACTORING SECTION
!  ---------------------------------------------------------------------

      call boundaries%fix_diagnostic_vars(tra_DIA, tra_DIA_2d)

! ----------------------------------------------------------------------
!  END BC_REFACTORING SECTION
!  ---------------------------------------------------------------------

#endif

                BIOparttime =  MPI_WTIME() -BIOparttime
                BIOtottime  = BIOtottime  + BIOparttime
               
      END SUBROUTINE trcbio
