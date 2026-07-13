SUBROUTINE trcbio_fabm(datestring)
#ifdef key_trc_fabm
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
      USE OPT_mem, ONLY: PAR, RMU
      USE BC_mem
      USE calendar
      USE mpi

!!!   FABM IMPLEMENTATION
      USE fabm


! ----------------------------------------------------------------------
!  BEGIN BC_REFACTORING SECTION
!  ---------------------------------------------------------------------

      use bc_set_mod

! ----------------------------------------------------------------------
!  END BC_REFACTORING SECTION
!  ---------------------------------------------------------------------
      
      IMPLICIT NONE

      character(LEN=17), INTENT(IN) ::  datestring
!!!----------------------------------------------------------------------
!!! local declarations
!!! ==================
      integer :: jk,jj,ji,jn
      INTEGER  :: year, month, day
      double precision :: sec
!     double precision :: yearday
!!!----------------------------------------------------------------------


call read_date_string(datestring, year, month, day, sec)


!   | --------------- |
!   | FABM MODEL CALL |
!   | --------------- |

BIOparttime = MPI_WTIME()

! Prepare all fields FABM needs to compute source terms (e.g., light)
yearday = DAY_OF_THE_YEAR(datestring)
!call model_fabm%link_scalar(fabm_standard_variables%number_of_days_since_start_of_the_year, yearday)
call model_fabm%prepare_inputs(SEC_FROM_START(datestring),year,month,day,sec)


! In the loops below, dy and w are local to the j,k point being processed.
! They would therefore need to be processed further within the loop to be included in
! the host's advection-diffusion-reaction treatment.
! Alternatively, the could be declared as (4D) global variables that are built up
! in the loop and processed after.

DO ji=1,jpi
   DO jj=1,jpj
      ! Retrieve tracer source terms (tracer units s-1).
      ! Array dy(1:nx, 1:size(model%interior_state_variables)) is assumed to be allocated.
      if (tmask(1,jj,ji) == 0) CYCLE
      call model_fabm%get_interior_sources(1, jpk, jj, ji, tra(:,jj,ji,:))
      ! Retrieve vertical velocities (sinking, floating, active movement) in m s-1.
      ! Array w(1:nx,1:size(model%interior_state_variables)) is assumed to be allocated.
      call model_fabm%get_vertical_movement(1, jpk, jj, ji, ogstm_sedipi(:,jj,ji,:))
!  Fluxes and source terms at the water surface
      flux_sf = 0.0D0
      sms_sf  = 0.0D0
      call model_fabm%get_surface_sources(jj, ji, flux_sf, sms_sf)

      tra(1,jj,ji,:) = tra(1,jj,ji,:) + flux_sf(:) / e3t(1,jj,ji) ! from flux to change in concentration (/ e3t)
   END DO
END DO

! Compute any remaining diagnostics
call model_fabm%finalize_outputs()


DO  jn=1,size(model_fabm%interior_diagnostic_variables)
      IF (model_fabm%interior_diagnostic_variables(jn)%save) THEN

         tra_DIA(jn)%data = model_fabm%get_interior_diagnostic_data(jn)
      END IF
END DO


DO jn=1,size(model_fabm%horizontal_diagnostic_variables)
      IF (model_fabm%horizontal_diagnostic_variables(jn)%save) THEN

         tra_DIA_2d(jn)%data = model_fabm%get_horizontal_diagnostic_data(jn)
      END IF
END DO

! ----------------------------------------------------------------------
!  BEGIN BC_REFACTORING SECTION
!  ---------------------------------------------------------------------

!       call boundaries%fix_diagnostic_vars(tra_DIA, tra_DIA_2d)

! ----------------------------------------------------------------------
!  END BC_REFACTORING SECTION
!  ---------------------------------------------------------------------


                BIOparttime =  MPI_WTIME() -BIOparttime
                BIOtottime  = BIOtottime  + BIOparttime
#endif         
      END SUBROUTINE trcbio_fabm
