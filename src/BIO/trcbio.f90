
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
      USE BC_mem
      USE mpi
      use mem, only: D3STATE

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
!      double precision,dimension(jpk,jptra) :: a
      double precision,dimension(4,jpk) :: c
      double precision,dimension(jptra_dia,jpk) :: d
      double precision,dimension(jpk,11) :: er
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
!      a        = 1.0
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

         er(jk,11) = correct_fact * ( gdept(jpk)-gdept(jk) ) /gdept(jpk)
      enddo
#endif


      DO ji=1,jpi
      DO jj=1,jpj
      if (bfmmask(1,jj,ji) == 0) CYCLE
      bottom = mbathy(jj,ji)

                          DO jtr=1, jtrmax

                             D3STATE(1:bottom, jtr) = trn(1:bottom,jj,ji,jtr) ! current biogeochemical concentrations

                          END DO

! Environmental regulating factors (er,:)

                          er(1:bottom,1)  = tn (1:bottom,jj,ji)! Temperature (Celsius)
                          er(1:bottom,2)  = sn (1:bottom,jj,ji)  ! Salinity PSU
                          er(1:bottom,3)  = rho(1:bottom,jj,ji)        ! Density Kg/m3
                          er(1       ,4)  = ice                  ! from 0 to 1 adimensional
                          er(1       ,5)  = ogstm_co2(jj,ji)     ! CO2 Mixing Ratios (ppm)  390
                          do jk=1, bottom
                          er(jk,6) = instant_par(COMMON_DATEstring,xpar(jk,jj,ji))  ! PAR umoles/m2/s | Watt to umoles photons W2E=1./0.217
                          enddo
                          !if (is_night(COMMON_DATEstring))  then
                          !    er(1:bottom,6)  = 0.001      
                          !else
                          !    er(1:bottom,6)  = 2.0*xpar(1:bottom,jj,ji)
                          !endif
                          !write(*,*) 'XPAR',  er(1,6)

                          er(1       ,7)  = DAY_LENGTH(jj,ji)    ! fotoperiod expressed in hours
                          er(1:bottom,8)  = e3t(1:bottom,jj,ji)        ! depth in meters of the given cell
                          er(1       ,9)  = vatm(jj,ji)                ! wind speed (m/s)
                          er(1:bottom,10) = ogstm_PH(1:bottom,jj,ji)   ! 8.1

#ifndef gdept1d
                         do jk=1,bottom
                             correct_fact= 1.0D0
                             if ( (gdept(jk,jj,ji) .GT. 1000.0D0 ) .AND.  (gdept(jk,jj,ji) .LT. 2000.0D0)) then
                                 correct_fact= 0.25D0
                             endif

                             if (gdept(jk,jj,ji) .GE. 2000.0D0 ) then
                                 correct_fact= 0.0D0
                             endif

                             er(jk,11) = correct_fact * ( gdept(jpk,jj,ji)-gdept(jk,jj,ji) ) /gdept(jpk,jj,ji)
                         enddo
#endif
                          call BFM1D_Input_EcologyDynamics(bottom,er)

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
