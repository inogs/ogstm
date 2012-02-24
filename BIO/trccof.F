C $Id: trccof.F,v 1.2 2009-09-11 09:20:56 cvsogs01 Exp $

      SUBROUTINE trccof 
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE trccof
CCC                     ******************
CCC
CCC  Purpose :
CCC  --------
CCC     Compute the arrays used for newtonian damping on passive tracer fields
CCC     ( if "key_trc_dmp" key activated )
CCC
CC   Method :
CC   -------
CC      Arrays defining the damping are computed for each grid
CC      point for each passive tracer (restotr)
CC      Damping depends on distance to coast, depth and latitude
CC
CC   Input :
CC   ------
CC      common
CC            /cotdmp/          :  passive tracer damping
CC
CC   Output :
CC   -------
CC      argument                : no
CC      common
CC            /cotdmp/          : passive tracer damping
CC
CC   External :
CC   --------- 
CC
CC   Modifications :
CC   --------------
CC      original  : 99-10 (H. Loukos) dtacof.F adapted for passive tracers
CC----------------------------------------------------------------------
CC parameters and commons
CC ======================

       USE myalloc
       USE myalloc_mpp
       USE BC_mem

       IMPLICIT NONE

CC----------------------------------------------------------------------
CC local declarations
CC ==================
C
#if defined key_trc_dmp && defined key_passivetrc
       INTEGER jn, jk, jj, ji, tra_idx

C
CC----------------------------------------------------------------------
CC statement functions
CC ===================

#include "stafun.h"

C
C===================================================
C===================================================
C
      DO jn=1,jptra
        DO jk=1,jpk-1
          DO jj=1,jpj
            DO ji=1,jpi
              restotr(ji,jj,jk,jn)=0.
            END DO
          END DO
        END DO
      END DO


      DO jn=1,jn_gib
        tra_idx=tra_matrix_gib(jn)
        DO jk=1,jpk
          DO jj=1,jpj
            DO ji=1,jpi
              restotr(ji,jj,jk,tra_idx)=resto(ji,jj,jk)
            END DO
          END DO
        END DO
      END DO
#  else
C
C     No evaluation of restoring coefficients 
C
#endif
C
C
      RETURN
      END
