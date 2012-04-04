      SUBROUTINE div ( ktask, kt )
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE div
CCC                     ***************
CCC
CCC  Purpose :
CCC  --------
CCC	compute the now horizontal divergence of the velocity field.
CCC
CC   Method :
CC   -------
CC    The now divergence is given by :
CC       * s-coordinate ('key_s_coord' defined)
CC         hdivn = 1/(e1t*e2t*e3t) ( di[e2u*e3u un] + dj[e1v*e3v vn] )
CC       * z-coordinate (default key)
CC         hdivn = 1/(e1t*e2t) [ di(e2u  un) + dj(e1v  vn) ]
CC
CC      Apply lateral boundary conditions on hdivn through a call to
CC      routine mpplnk ('key_mpp' defined) or lbc.
CC
CC      macro-tasked on horizontal slab (jk-loop :  1, jpk-1)
CC
CC   Input :
CC   ------
CC      argument
CC              ktask           : task identificator
CC              kt              : time step
CC      common
CC            /comcoh/       : scale factors
CC            /comtsk/       : multitasking
CC            /comnow/        : present fields (now)
CC
CC   Output :
CC   -------
CC      common
CC          /comnow/ hdivn    : now horizontal divergence
CC
CC   External : mpplnk or lbc

      USE myalloc
      USE stafun
      IMPLICIT NONE
CC    local declarations

      INTEGER ktask, kt
      INTEGER ji, jj, jk

      REAL(8) zbt
      REAL(8) zwu(jpi,jpj), zwv(jpi,jpj)
CC----------------------------------------------------------------------
CC statement functions
CC ===================

!                                          #include "stafun.h"

C Horizontal slab
C ===============

      DO jk = ktask, jpkm1


C 1. Horizontal fluxes
C --------------------
C
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
#if defined key_s_coord 
            zwu(ji,jj) = e2u(ji,jj) * fse3u(ji,jj,jk) * un(ji,jj,jk)
            zwv(ji,jj) = e1v(ji,jj) * fse3v(ji,jj,jk) * vn(ji,jj,jk)
#  else
            zwu(ji,jj) = e2u(ji,jj) * un(ji,jj,jk)
            zwv(ji,jj) = e1v(ji,jj) * vn(ji,jj,jk)
#endif
          END DO
        END DO
C
C
C 2. horizontal divergence
C ------------------------
C
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zbt = e1t(ji,jj) * e2t(ji,jj)
#if defined key_s_coord
     $                                    * fse3t(ji,jj,jk)
#endif
            hdivn(ji,jj,jk) = (  zwu(ji,jj) - zwu(ji-1,jj  )
     $                         + zwv(ji,jj) - zwv(ji  ,jj-1)  ) / zbt
          END DO  
        END DO  
C
C
C End of slab
C ===========
C
      END DO
C
C
C 3. Lateral boundary conditions on hdivn
C ---------------------------------=======
#ifdef key_mpp
C
C ... Mpp : export boundary values to neighboring processors
      CALL mpplnk( hdivn, 1, 1 )
#  else
C
C ... mono or macro-tasking: T-point, 3D array, jk-slab
      CALL lbc( hdivn, 1, 1, 1, ktask, jpkm1, 1 )
#endif
C
C
      RETURN
      END
