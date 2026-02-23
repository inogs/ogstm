
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

      double precision,dimension(jpk,jptra) :: dy
      double precision,dimension(4,jpk) :: c
!     double precision,dimension(jptra_dia,jpk) :: d
      double precision,dimension(jpk,16) :: er
      double precision,dimension(jptra_dia_2d) :: d2


      integer :: jk,jj,ji,jb,jn
      integer :: ivar
      integer :: jtr,jtrmax,tra_idx
      integer :: bottom
      double precision :: correct_fact


!!!----------------------------------------------------------------------
!!! ===================


!   | --------------|
!   | FABM MODEL CALL|
!   | --------------|

!       BIOparttime = MPI_WTIME()


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
      call model_fabm%get_interior_sources(1, jk, jj, ji, tra(:,jj,ji,:))
      ! Retrieve vertical velocities (sinking, floating, active movement) in m s-1.
      ! Array w(1:nx,1:size(model%interior_state_variables)) is assumed to be allocated.
!      call model_fabm%get_vertical_movement(1, jpk, jj, ji, w)
   end do
end do

! Compute any remaining diagnostics
call model_fabm%finalize_outputs()
!         surf_mask(:) = 0.
!         surf_mask(1) = 1.
! -------------------------------------------------

          ! tra_idx = tra_matrix_gib(1)
!         jtrmax=jptra

! ---------------- Fuori dai punti BFM

!     ogstm_sediPI=0.
!     tra_DIA    = 0.
!     tra_DIA_2d = 0. ! da sistemare


!    Initialization
!     a        = 1.0
!     er       = 1.0
!     er(:,10) = 8.1


!     DO ji=1,jpi
!     DO jj=1,jpj
!     if (bfmmask(1,jj,ji) == 0) CYCLE
!     bottom = mbathy(jj,ji)

!                         DO jtr=1, jtrmax

!                            a(1:bottom, jtr) = trn(1:bottom,jj,ji,jtr) ! current biogeochemical concentrations

!                         END DO

! Environmental regulating factors (er,:)


!                         call BFM1D_Input_EcologyDynamics(bottom,a,jtrmax,er)

!                        call BFM1D_reset()

!                        call EcologyDynamics()

!                        call BFM1D_Output_EcologyDynamics(b, c, d, d2)

!                         DO jtr=1, jtrmax
!                            tra(1:bottom,jj,ji,jtr) =tra(1:bottom,jj,ji,jtr) +b(jtr,1:bottom) ! trend
!                         END DO

!                         DO jtr=1,4
!                            ogstm_sediPI(1:bottom,jj,ji,jtr) = c(jtr,1:bottom)      ! BFM output of sedimentation speed (m/d)
!                         END DO


!                         DO jk = 1,bottom
!                         DO jtr=1,jptra_dia
!                            tra_DIA(jtr, jk ,jj,ji) = d(jtr,jk) ! diagnostic
!                         END DO
!                         ENDDO

!                        tra_DIA_2d(:,jj,ji) = d2(:) ! diagnostic

!                        ogstm_PH(1:bottom,jj,ji) = d(pppH,1:bottom) ! Follows solver guess, put 8.0 if pppH is not defined

!     END DO
!     END DO

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
