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

c#include "parameter.h"
c#include "common.h"
       USE myalloc
       USE myalloc_mpp

CC+CC  Implicit typing is never allowed
        IMPLICIT NONE
CC+CC  Implicit typing is never allowed

CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER kkflag
C
#if defined key_trc_dmp && defined key_passivetrc
       INTEGER jn, jk, jj, ji, jpp, irest
       INTEGER jn_nudg,tra_matrix(51)
       REAL(8) rs_ts(13), fact, fact1, tmaxt

C
CC----------------------------------------------------------------------
CC statement functions
CC ===================

#include "stafun.h"

CCC---------------------------------------------------------------------
CCC  OPA8, LODYC (15/11/96)
CCC---------------------------------------------------------------------
      jn_nudg = 4
C
      tra_matrix(1) = ppO2o
      tra_matrix(2) = ppN1p
      tra_matrix(3) = ppN3n
      tra_matrix(4) = ppN5s
C
C

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


      DO jn=1,jn_nudg
        DO jk=1,jpk
          DO jj=1,jpj
            DO ji=1,jpi
              restotr(ji,jj,jk,tra_matrix(jn))=resto(ji,jj,jk)
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
