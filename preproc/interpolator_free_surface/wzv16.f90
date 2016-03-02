      SUBROUTINE wzv16 ( )
!CC---------------------------------------------------------------------
!CC
!CC                       ROUTINE wzv
!CC                     ***************
!CC
!CC  Purpose :
!CC  ---------
!CC	Compute the now vertical velocity after the array swap.
!CC
!C   Method :
!C   -------
!C	Using the incompressibility hypothesis, the vertical velocity
!C	is computed by integrating the horizontal divergence from the
!C	bottom to the surface.
!C	The boundary conditions are w=0 at the bottom (no flux) and
!C	w=0 at the sea surface (rigid lid).
!C
!C	macro-tasked on vertical slab (jj-loop)
!C
!C   Input :
!C   ------
!C      argument
!C              ktask       : task identificator
!C              kt          : time step
!C      common
!C              /comcoh/    : scale factors
!C		/comtsk/    : multitasking
!C              /comnow/    : present fields (now)
!C
!C   Output :
!C   -------
!C      common
!C	       /comnow/ wn  : now vertical velocity 
!C
!C   Modifications :
!C   --------------
!C      original : 90-10 (cl-gm)
!C      addition : 91-11 (G. Madec)
!C      addition : 96-01 (G. Madec) statement function for e3
!C----------------------------------------------------------------------
!C parameters and commons
!C ======================
       USE modulo16
!C+CC  Implicit typing is never allowed
        IMPLICIT NONE
!C+CC  Implicit typing is never allowed

!C----------------------------------------------------------------------
!C local declarations
!C ==================
      INTEGER ktask, kt
      INTEGER ji16, jj16, jk
      INTEGER jpi16m1,jpj16m1,jpk16m1

!C----------------------------------------------------------------------
!C statement functions
!C ===================
!CC---------------------------------------------------------------------
!CC  OPA8, LODYC (1997)
!CC---------------------------------------------------------------------
!
!
! Vertical slab
! =============
!
      jpi16m1=jpi16-1
      jpj16m1=jpj16-1
      jpk16m1=jpk16-1
      DO 1000 jj16 = 1, jpj16
!
!
! 1. Surface and bottom boundary condition: w=0 (rigid lid and no flux)
! ----------------------------------------
!
        DO ji16 = 1, jpi16
          wn16(ji16,jj16, 1 ) = 0.e0
          wn16(ji16,jj16,jpk16) = 0.e0
        END DO  
!
!
! 2. Computation from the bottom
! ------------------------------
!
         DO jk = jpk16m1, 1, -1
          DO ji16 = 1, jpi16
            wn16(ji16,jj16,jk) = wn16(ji16,jj16,jk+1)  &
                        - e3t16(ji16,jj16,jk)*hdivn16(ji16,jj16,jk)
          END DO 
        END DO
!
!
! End of slab
! ===========
!
 1000 CONTINUE
!
!
      RETURN
      END
