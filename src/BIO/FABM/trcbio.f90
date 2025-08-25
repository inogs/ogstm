
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


!!!   FABM IMPLEMENTATION
      USE fabm
!  ---------------------------------------------------------------------

! ----------------------------------------------------------------------
!  BEGIN BC_REFACTORING SECTION
!  ---------------------------------------------------------------------

      use bc_set_mod

! ----------------------------------------------------------------------
!  END BC_REFACTORING SECTION
!  ---------------------------------------------------------------------
      
#ifdef ExecEnsParams
      use Ens_Mem, &
          only: UseParams
      use Ens_Params, &
          only: Ens_SetParams_trcbio
#endif
      
      IMPLICIT NONE


!!!----------------------------------------------------------------------
!!! local declarations
!!! ==================

      double precision,dimension(jptra,jpk) :: b
      double precision,dimension(jpk,jptra) :: a
      double precision,dimension(4,jpk) :: c
      double precision,dimension(jptra_dia,jpk) :: d
      double precision,dimension(jpk,11) :: er
      double precision,dimension(jptra_dia_2d) :: d2


      integer :: jk,jj,ji,jb,jn
      integer :: jtr,jtrmax,tra_idx
      integer :: bottom
      double precision :: correct_fact

!!! FABM IMPLEMENTATION
      class (type_fabm_model), pointer :: model

      ! Initialize (reads FABM configuration from fabm.yaml)
      ! After this the number of biogeochemical variables is fixed.
      ! (access variable metadata in model%interior_state_variables, model%interior_diagnostic_variables)
      model => fabm_create_model()
      
!  ---------------------------------------------------------------------


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

         er(jk,11) = correct_fact * ( gdept(jpk)-gdept(jk) ) /gdept(jpk)
      enddo
#endif


      DO ji=1,jpi
      DO jj=1,jpj
      if (bfmmask(1,jj,ji) == 0) CYCLE
      bottom = mbathy(jj,ji)
      
!!! FABM IMPLEMENTATION
      ! Provide extents of the spatial domain (number of layers nz for a 1D column)
!     call model%set_domain(nz)  !maybe nz is jpk or bottom in ogstm
      call model%set_domain(jpk) 

! At this point (after the call to fabm_create_model), memory should be
! allocated to hold the values of all size(model%interior_state_variables) state variables.
! Where this memory resides and how it is laid out is typically host-specific.
! Below, we assume all state variable values are combined in an array interior_state with
! shape nz,size(model%interior_state_variables).

      ! Point FABM to your state variable data
!     real, dimension(:,:), target, allocatable :: interior_state
!     allocate(interior_state        (nz, size(model%interior_state_variables))
      ! interior_state is "a" in the original ogstm-bfm coupling
      ! size(model%interior_state_variables) is "jptra" in the original ogstm-bfm coupling
!     do ivar = 1, size(model%interior_state_variables)
      DO ivar = 1, jptra
         call model%link_interior_state_data(ivar, a(ivar,:))  ! current biogeochemical concentrations
      END DO
!  ---------------------------------------------------------------------

#ifdef ExecEnsParams
      if (UseParams) call Ens_SetParams_trcbio(jj, ji)
#endif


!                         DO jtr=1, jtrmax

!                            a(1:bottom, jtr) = trn(1:bottom,jj,ji,jtr) ! current biogeochemical concentrations

!                         END DO

! Environmental regulating factors (er,:)

!!! FABM IMPLEMENTATION

! Point FABM to environmental data, here shown for temperature
! Array temp(1:nz) is assumed to be allocated.
! Do this for all variables on FABM's standard variable list that the model can provide.
! For this list, visit https://fabm.net/standard_variables
! call model%link_interior_data(fabm_standard_variables%temperature, temp) example

! if in a loop better another approach, see wiki

call model%link_interior_data(fabm_standard_variables%temperature, tn(1:bottom,jj,ji))        ! Temperature (Celsius)
call model%link_interior_data(fabm_standard_variables%practical_salinity, sn(1:bottom,jj,ji)) ! salinity PSU
call model%link_interior_data(fabm_standard_variables%density, rho(1:bottom,jj,ji))           ! Density Kg/m3
call model%link_horizontal_data(fabm_standard_variables%ice_area_fraction, ice)                   ! from 0 to 1 adimensional
call model%link_horizontal_data(fabm_standard_variables%mole_fraction_of_carbon_dioxide_in_air, ogstm_co2(jj,ji))                   ! CO2 Mixing Ratios (ppm)  390
call model%link_interior_data(fabm_standard_variables%depth, e3t(1:bottom,jj,ji))             ! depth in meters of the given cell
call model%link_horizontal_data(fabm_standard_variables%wind_speed, vatm(jj,ji))              ! depth in meters of the given cell
call model%link_interior_data(fabm_standard_variables%ph_reported_on_total_scale, ogstm_PH(1:bottom,jj,ji))             ! 8.1
!how to do it in fabm?
!do jk=1, bottom
!    er(jk,6) = instant_par(COMMON_DATEstring,xpar(jk,jj,ji))  ! PAR umoles/m2/s | Watt to umoles photons W2E=1./0.217
!enddo
!er(1       ,7)  = DAY_LENGTH(jj,ji)    ! fotoperiod expressed in hours

!  ---------------------------------------------------------------------

!                         er(1:bottom,1)  = tn (1:bottom,jj,ji)! Temperature (Celsius)
!                         er(1:bottom,2)  = sn (1:bottom,jj,ji)  ! Salinity PSU
!                         er(1:bottom,3)  = rho(1:bottom,jj,ji)        ! Density Kg/m3
!                         er(1       ,4)  = ice                  ! from 0 to 1 adimensional
!                         er(1       ,5)  = ogstm_co2(jj,ji)     ! CO2 Mixing Ratios (ppm)  390
!                         do jk=1, bottom
!                         er(jk,6) = instant_par(COMMON_DATEstring,xpar(jk,jj,ji))  ! PAR umoles/m2/s | Watt to umoles photons W2E=1./0.217
!                         enddo
!                         !if (is_night(COMMON_DATEstring))  then
!                         !    er(1:bottom,6)  = 0.001      
!                         !else
!                         !    er(1:bottom,6)  = 2.0*xpar(1:bottom,jj,ji)
!                         !endif
!                         !write(*,*) 'XPAR',  er(1,6)

!                         er(1       ,7)  = DAY_LENGTH(jj,ji)    ! fotoperiod expressed in hours
!                         er(1:bottom,8)  = e3t(1:bottom,jj,ji)        ! depth in meters of the given cell
!                         er(1       ,9)  = vatm(jj,ji)                ! wind speed (m/s)
!                         er(1:bottom,10) = ogstm_PH(1:bottom,jj,ji)   ! 8.1

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

!!! FABM IMPLEMENTATION
                         ! Complete initialization and check whether FABM has all dependencies fulfilled
                         ! (i.e., whether all required calls to model%link_*_data have been made)
                         call model%start()
                         
                         ! Initialize the tracers
                         ! This sets the values of arrays sent to model%link_interior_state_data,
                         ! in this case those contained in interior_state.
                         call model%initialize_interior_state(1, jpk)
                         
                         ! At this point, initialization is complete.
                         ! Routines below would typically be called every time step.
                         
                         ! Prepare all fields FABM needs to compute source terms (e.g., light)
                         call model%prepare_inputs()
                         
                         ! Retrieve tracer source terms (tracer units s-1).
                         ! Array dy(1:nz,1:size(model%interior_state_variables)) is assumed to be allocated.
                         !dy = 0
                         !call model%get_interior_sources(1, nz, dy)
!                        real(rk) :: b(jpk,size(model%interior_state_variables))
                         real(rk) :: b(jtrmax,jpk)
                         b = 0
                         call model%get_interior_sources(1, jpk, b)
                         
                         ! Compute any remaining diagnostics
                         call model%finalize_outputs()


                         ! Retrieve vertical velocities (sinking, floating, active movement) in m s-1.
                         ! Array w(1:nz, 1:size(model%interior_state_variables)) is assumed to be allocated.

                         !call model%get_vertical_movement(1, nz, w) not needed in bfm?
                         ! Here you would time-integrate the advection-diffusion-reaction equations
                         ! of all tracers, combining the transport terms with the biogeochemical source
                         ! terms dy and vertical velocities w. This should result in an updated interior_state.

!  ---------------------------------------------------------------------


!                         call BFM1D_Input_EcologyDynamics(bottom,a,jtrmax,er)

!                        call BFM1D_reset()

!                        call EcologyDynamics()

!                        call BFM1D_Output_EcologyDynamics(b, c, d, d2)

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
