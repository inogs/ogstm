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
       IMPLICIT NONE

#include <mpif.h>

      trcstpparttime = MPI_WTIME() !  F79 cronometer-start

! this part is parallelized tracer  by tracer
      CALL trcadv ! tracers: advection


#    if defined key_trc_dmp
        CALL trcdmp ! tracers: damping for passive tracers
#    else
!       no passive tracer damping
#    endif

! tracers: horizontal diffusion IF lhdf is true
! -----------------------------

      IF (lhdf)            CALL trchdf
      IF ( ntrbbl .EQ. 1 ) CALL trcbbl


! tracers: sink and source (must be  parallelized on vertical slab)
! ------------------------
      CALL trcsbc ! surface cell processes
      CALL trcsms
      CALL trczdf ! tracers: vertical diffusion
      CALL snutel
      CALL checkValues
      CALL trcnxt ! tracers: fields at next time step

      trcstpparttime = MPI_WTIME() - trcstpparttime ! F79 cronometer-stop
      trcstptottime = trcstptottime + trcstpparttime

      END SUBROUTINE trcstp

