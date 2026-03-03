SUBROUTINE trcbio_fabm
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


!!!----------------------------------------------------------------------
!!! local declarations
!!! ==================
      integer :: jk,jj,ji,jn
      double precision,dimension(jpk,jptra) :: dy
!!!----------------------------------------------------------------------



!   | --------------- |
!   | FABM MODEL CALL |
!   | --------------- |

       BIOparttime = MPI_WTIME()


! Prepare all fields FABM needs to compute source terms (e.g., light)
call model_fabm%prepare_inputs()

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
      dy = 0
      call model_fabm%get_interior_sources(1, jpk, jj, ji, tra(:,jj,ji,:))
      ! Retrieve vertical velocities (sinking, floating, active movement) in m s-1.
      ! Array w(1:nx,1:size(model%interior_state_variables)) is assumed to be allocated.
       call model_fabm%get_vertical_movement(1, jpk, jj, ji, ogstm_sedipi(:,jj,ji,:))
   END DO
END DO

! Compute any remaining diagnostics
call model_fabm%finalize_outputs()


DO  jn=1,size(model_fabm%interior_diagnostic_variables)
      IF (model_fabm%interior_diagnostic_variables(jn)%save) THEN
!      write(*,*) 'jn=', jn
!      write(*,*) 'Shape of get_interior_diagnostic_data(jn):', shape(model_fabm%get_interior_diagnostic_data(jn))
!      write(*,*) 'Shape of tra_DIA(jn,:,:,:):', shape(tra_DIA(jn,:,:,:))
         tra_DIA(jn, : , : ,:) = model_fabm%get_interior_diagnostic_data(jn)
      ELSE
         tra_DIA(jn, : , :, :) = 0.
      END IF
END DO


DO jn=1,size(model_fabm%horizontal_diagnostic_variables)
      IF (model_fabm%horizontal_diagnostic_variables(jn)%save) THEN
         tra_DIA_2d(jn, :, :) = model_fabm%get_horizontal_diagnostic_data(jn)
      ELSE
         tra_DIA_2d(jn, :, :) = 0.
      END IF
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
#endif         
      END SUBROUTINE trcbio_fabm
