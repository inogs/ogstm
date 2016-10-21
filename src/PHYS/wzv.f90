      SUBROUTINE wzv ()
!---------------------------------------------------------------------
!
!                       ROUTINE wzv
!                     ***************
!
!  Purpose :
!  ---------
!	Compute the now vertical velocity after the array swap.
!
!   Method :
!   -------
!	Using the incompressibility hypothesis, the vertical velocity
!	is computed by integrating the horizontal divergence from the
!	bottom to the surface.
!	The boundary conditions are w=0 at the bottom (no flux) and
!	w=0 at the sea surface (rigid lid).
!
!	macro-tasked on vertical slab (jj-loop)
!
!   Input :
!   ------
!      argument
!              ktask       : task identificator
!              kt          : time step
!      common
!              /comcoh/    : scale factors
!		/comtsk/    : multitasking
!              /comnow/    : present fields (now)
!
!   Output :

! parameters and commons
! ======================

       USE myalloc
       IMPLICIT NONE

!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER ji, jj, jk

! Vertical slab
! =============

      DO jj = 1, jpj

! 1. Surface and bottom boundary condition: w=0 (rigid lid and no flux)
! ----------------------------------------

        DO ji = 1, jpi
          wdta(jj,ji, 1 ,2) = 0.e0
          wdta(jpk,jj,ji,2) = 0.e0
        END DO  


! 2. Computation from the bottom
! ------------------------------

        DO jk = jpkm1, 1, -1
          DO ji = 1, jpi
            wdta(jk,jj,ji,2) = wdta(jk+1,jj,ji,2) - e3t(jk,jj,ji)*hdivn(jk,jj,ji)
!           wn(jk,jj,ji) = wn(jk,jj,ji+1) - e3t(jk,jj,ji)*hdivn(jk,jj,ji)
          END DO 
        END DO

      END DO

      END SUBROUTINE wzv
