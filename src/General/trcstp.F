      SUBROUTINE trcstp
!---------------------------------------------------------------------
!
!                       ROUTINE trcstp
!                     *****************
!
!  PURPOSE :
!  ---------
!	time loop of ogstm for passive tracer
!
!   METHOD :
!   -------
!      compute the well/spring evolution
!
!      compute the time evolution of tracers concentration
!         with advection
!         with horizontal diffusion
!         with surface boundary condition
!         with IMPLICIT vertical diffusion


       USE myalloc
       USE myalloc_mpp
       IMPLICIT NONE

      trcstpparttime = MPI_WTIME() ! cronometer-start

      CALL trcadv ! tracers: advection

#    if defined key_trc_dmp
      CALL trcdmp ! tracers: damping for passive tracers
#    endif

! tracers: horizontal diffusion IF namelist flags are activated
! -----------------------------

      IF (lhdf)   CALL trchdf


! tracers: sink and source (must be  parallelized on vertical slab)
! ------------------------
      CALL trcsbc ! surface cell processes
      CALL trcsms
      CALL trczdf ! tracers: vertical diffusion
      CALL snutel
      CALL hard_tissue_pump
      CALL checkValues
      CALL trcnxt ! tracers: fields at next time step

      trcstpparttime = MPI_WTIME() - trcstpparttime ! cronometer-stop
      trcstptottime = trcstptottime + trcstpparttime

      END SUBROUTINE trcstp

