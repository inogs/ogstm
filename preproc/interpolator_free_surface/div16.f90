      SUBROUTINE div16 (jk)
!CC---------------------------------------------------------------------
!CC
!CC                       ROUTINE div
!CC                     ***************
!CC
!CC  Purpose :
!CC  --------
!CC	compute the now horizontal divergence of the velocity field.
!CC
!C   Method :
!C   -------
!C	The now divergence is given by :
!C       * s-coordinate ('key_s_coord' defined)
!C         hdivn = 1/(e1t*e2t*e3t) ( di[e2u*e3u un] + dj[e1v*e3v vn] )
!C       * z-coordinate (default key)
!C         hdivn = 1/(e1t*e2t) [ di(e2u  un) + dj(e1v  vn) ]
!C
!C      Apply lateral boundary conditions on hdivn through a call to
!C      routine mpplnk ('key_mpp' defined) or lbc.
!C
!C      macro-tasked on horizontal slab (jk-loop :  1, jpk-1)
!C
!C   Input :
!C   ------
!C      argument
!C              ktask           : task identificator
!C              kt              : time step
!C      common
!C            /comcoh/   	: scale factors
!C            /comtsk/   	: multitasking
!C            /comnow/		: present fields (now)
!C
!C   Output :
!C   -------
!C      common
!C	      /comnow/ hdivn	: now horizontal divergence
!C
!C   External : mpplnk or lbc
!C   ---------
!C
!C   Modifications :
!C   --------------
!C      original :  87-06 (P. Andrich, D. L Hostis)
!C      additions : 91-11 (G. Madec)
!C                : 93-03 (M. Guyon) symetrical conditions
!C                : 96-01 (G. Madec) s-coordinates
!C                : 97-06 (G. Madec) lateral boundary cond., lbc
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
      INTEGER ji16, jj16, jk
      INTEGER jpi16m1,jpj16m1,jpk16m1
!
      REAL(4) zbt
      REAL(4) zwu16(jpi16,jpj16), zwv16(jpi16,jpj16)
!C----------------------------------------------------------------------
!C statement functions
!C ===================
!CC---------------------------------------------------------------------
!CC  OPA8, LODYC (1997)
!CC---------------------------------------------------------------------
!
! Horizontal slab
! ===============
!
      jpi16m1=jpi16-1
      jpj16m1=jpj16-1
      jpk16m1=jpk16-1
      
!C-CC      DO 1000 jk = 1, jpk16m1 
!
!
! 1. Horizontal fluxes
! --------------------
!
        DO jj16 = 1, jpj16m1
          DO ji16 = 1, jpi16m1
!C-CC            zwu16(ji16,jj16) = e2u16(ji16,jj16) * fse3u16(ji16,jj16,jk) * un16(ji16,jj16,jk)
!C-CC            zwv16(ji16,jj16) = e1v16(ji16,jj16) * fse3v16(ji16,jj16,jk) * vn16(ji16,jj16,jk)

            zwu16(ji16,jj16) = e2u16(ji16,jj16) * e3u16(ji16,jj16,jk) * un16(ji16,jj16,jk)
            zwv16(ji16,jj16) = e1v16(ji16,jj16) * e3v16(ji16,jj16,jk) * vn16(ji16,jj16,jk)

          END DO
        END DO
!
!
! 2. horizontal divergence
! ------------------------
!
        DO jj16 = 2, jpj16m1
          DO ji16 = 2, jpi16m1
            zbt = e1t16(ji16,jj16) * e2t16(ji16,jj16)* e3t16(ji16,jj16,jk)
            hdivn16(ji16,jj16,jk) = (  zwu16(ji16,jj16) - zwu16(ji16-1,jj16  )  &
                                + zwv16(ji16,jj16) - zwv16(ji16  ,jj16 -1)  ) / zbt
          END DO  
        END DO  
!
!
! End!of slab
! ===========
!
!C-CC 1000 CONTINUE
!
!
! 3. Lateral boundary conditions on hdivn
! ---------------------------------=======
!C-CC#ifdef key_mpp
!
! ... Mpp : export boundary values to neighboring processors
!C-CC      CALL mpplnk( hdivn, 1, 1 )
!C-CC#  else
!
! ... mono or macro-tasking: T-point, 3D array, jk-slab
!C-CC      CALL lbc( hdivn, 1, 1, 1, ktask, jpkm1, ncpu )
!C-CC#endif
!
!
      RETURN
      END
