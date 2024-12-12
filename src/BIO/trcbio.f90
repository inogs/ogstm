
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
      USE OPT_mem, ONLY: PAR, RMU,SWR_RT
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

      double precision,dimension(jptra,jpk) :: b
      double precision,dimension(jpk,jptra) :: a
      double precision,dimension(4,jpk) :: c
      double precision,dimension(jptra_dia,jpk) :: d
      double precision,dimension(jpk,18) :: er
      double precision,dimension(jptra_dia_2d) :: d2


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

                          er(1:bottom,1)  = tn (1:bottom,jj,ji)! Temperature (Celsius)
                          er(1:bottom,2)  = sn (1:bottom,jj,ji)  ! Salinity PSU
                          er(1:bottom,3)  = rho(1:bottom,jj,ji)        ! Density Kg/m3
                          er(1       ,4)  = ice                  ! from 0 to 1 adimensional
                          er(1       ,5)  = ogstm_co2(jj,ji)     ! CO2 Mixing Ratios (ppm)  390
                          er(1:bottom,6)  = PAR(1:bottom,jj,ji,1) ! PAR for diatoms
                          er(1:bottom,7)  = PAR(1:bottom,jj,ji,2) ! PAR for flagellates
                          er(1:bottom,8)  = PAR(1:bottom,jj,ji,3) ! PAR for pico phytoplankton
                          er(1:bottom,9)  = PAR(1:bottom,jj,ji,4) ! PAR for dinoflagellates
                          er(1:bottom,10) = PAR(1:bottom,jj,ji,5) ! total PAR for CDOM
                          er(1       ,11)  = DAY_LENGTH(jj,ji)    ! fotoperiod expressed in hours
                          er(1:bottom,12)  = e3t(1:bottom,jj,ji)        ! depth in meters of the given cell
                          er(1       ,13)  = vatm(jj,ji)                ! wind speed (m/s)
                          er(1:bottom,14) = ogstm_PH(1:bottom,jj,ji)   ! 8.1
                          er(1       ,15)  = RMU(jj,ji)                ! avg. cosine direct

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

                            er(1:bottom,17) = SWR_RT(1:bottom,jj,ji)  
                            er(1       ,18)  = atm_Hg0(jj,ji)         
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


                BIOparttime =  MPI_WTIME() -BIOparttime
                BIOtottime  = BIOtottime  + BIOparttime
               
      END SUBROUTINE trcbio
