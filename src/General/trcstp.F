      SUBROUTINE trcstp(kt)
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
!
!   INPUT :
!   -----
!      argument
!              kt              : time step
!   EXTERNAL :
!   --------
!      trc3ad, trcsbc, trczdf, trchdf, trcnxt, trcsms
!                        passive tracers interface

       USE myalloc
       IMPLICIT NONE

#include <mpif.h>


!     local declarations
      INTEGER kt
      INTEGER ktask


      ktask = 1
      trcstpparttime = MPI_WTIME() !  F79 cronometer-start


! this ROUTINE is called only every ndttrc time step
      IF (mod(kt,ndttrc) .NE. 0) RETURN


! this part is parallelized tracer  by tracer


      CALL trcadv(ktask,kt) ! tracers: advection


#    if defined key_trc_dmp
        CALL trcdmp ! tracers: damping for passive tracers
#    else
!       no passive tracer damping
#    endif

! tracers: horizontal diffusion IF lhdf is true
! -----------------------------

      IF (lhdf)            CALL trchdf(ktask,kt)
      IF ( ntrbbl .EQ. 1 ) CALL trcbbl(ktask,kt)


! tracers: sink and source (must be  parallelized on vertical slab)
! ------------------------
      CALL trcsbc ! surface cell processes
      CALL trcsms(ktask,kt)
      CALL trczdf(ktask,kt) ! tracers: vertical diffusion
      CALL snutel
      CALL checkValues
      CALL trcnxt(ktask,kt) ! tracers: fields at next time step
!      CALL trcrad           ! Correct artificial negative concentrations for isopycnal scheme
!
!
      trcstpparttime = MPI_WTIME() - trcstpparttime ! F79 cronometer-stop
      trcstptottime = trcstptottime + trcstpparttime

      END SUBROUTINE trcstp






!   ********* SEDIMENT MASS CONTROL *************************

! Controle masse dans le sediment
      SUBROUTINE SEDIMENT_MASS_CONTROL

      USE myalloc
      IMPLICIT NONE
      REAL(8) zmasspoc,zmasssil,zmasscal
      INTEGER ji,jj

      real(8) sedpoc(2,2), sedcal(2,2), sedsil(2,2) ! *** shoud be defined elsewhere !!!! ***


!      if (mod(kt,1095).eq.0) then

       zmasspoc = 0.
       zmasscal = 0.
       zmasssil = 0.

       DO jj = 2,jpj-1
          DO ji = 2,jpi-1
             zmasspoc = zmasspoc + 12/1e12*sedpoc(ji,jj)*tmask(ji,jj,1)*e1t(ji,jj)*e2t(ji,jj)
             zmasscal = zmasscal + 12/1e12*sedcal(ji,jj)*tmask(ji,jj,1)*e1t(ji,jj)*e2t(ji,jj)
             zmasssil = zmasssil + 1/1e9  *sedsil(ji,jj)*tmask(ji,jj,1)*e1t(ji,jj)*e2t(ji,jj)
          ENDDO
       ENDDO

       write(0,*) 'Total Mass for POC, CAL, SIL',zmasspoc,zmasscal,zmasssil

!       endif

      END SUBROUTINE SEDIMENT_MASS_CONTROL
