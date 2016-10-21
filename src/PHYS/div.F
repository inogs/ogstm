
      SUBROUTINE div ()
!!!---------------------------------------------------------------------
!!!
!!!                       ROUTINE div
!!!                     ***************
!!!
!!!  Purpose :
!!!  --------
!!!	compute the now horizontal divergence of the velocity field.
!!!
CC   Method :
CC   -------
CC	The now divergence is given by :
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
CC            /comcoh/   	: scale factors
CC            /comtsk/   	: multitasking
CC            /comnow/		: present fields (now)
CC
CC   Output :
CC   -------
CC      common
CC	      /comnow/ hdivn	: now horizontal divergence
CC
CC   External : mpplnk or lbc

      USE myalloc
      IMPLICIT NONE
CC    local declarations

      INTEGER ktask, kt
      INTEGER jk,jj,ji, ijt

      REAL(8) zbt
      REAL(8) zwu(jpj,jpi), zwv(jpj,jpi)
CC----------------------------------------------------------------------
CC statement functions
CC ===================

C Horizontal slab
C ===============

      DO jk = 1, jpkm1


C 1. Horizontal fluxes
C --------------------
C
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
#if defined key_s_coord 
            zwu(jj,ji) = e2u(jj,ji) * fse3u(jk,jj,ji) * udta(jk,jj,ji,2)
            zwv(jj,ji) = e1v(jj,ji) * fse3v(jk,jj,ji) * vdta(jk,jj,ji,2)
#  else
            zwu(jj,ji) = e2u(jj,ji) * e3u(jk,jj,ji) * udta(jk,jj,ji,2)
            zwv(jj,ji) = e1v(jj,ji) * e3v(jk,jj,ji) * vdta(jk,jj,ji,2)
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
            zbt = e1t(jj,ji) * e2t(jj,ji) * e3t(jk,jj,ji)
#if defined key_s_coord
     $                                    * fse3t(jk,jj,ji)
#endif
            hdivn(jk,jj,ji) = (  zwu(jj,ji) - zwu(jj,ji-1  )
     $                         + zwv(jj,ji) - zwv(jj-1  ,ji)  ) / zbt
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
!     CALL mpplnk_my( hdivn,1, 1, 1 )
#  else
C
C ... mono or macro-tasking: T-point, 3D array, jk-slab
!     CALL lbc( hdivn, 1, 1, 1, ktask, jpkm1, 1 )
#endif
C
C
      RETURN
      END
