C $Id: trcdtr.F,v 1.2 2009-09-11 09:20:56 cvsogs01 Exp $

      SUBROUTINE trcdtr
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE trcdtr
CCC                     ******************
CCC
CCC  PURPOSE :
CCC  ---------
CCC     computes or READ initial DATA for passive tracer
CCC
CC   METHOD :
CC   -------
CC
CC   INPUT :			  no
CC   -----
CC      COMMON
CC            /comdom/          : domain PARAMETER
CC            /comcoo/          : orthogonal curvilinear coordinates
CC                                and scale factors
CC            /comask/          : masks, bathymetry
CC
CC   OUTPUT :
CC   ------
CC      argument                : no
CC      COMMON
CC            /cottrc/          : passive tracer field now and before
CC
CC   WORKSPACE :		  no
CC   ---------
CC
CC   MODIFICATIONS:
CC   --------------
CC      original  : 96-11
CC      additions : 99-9
CC                : 00-12 (O. Aumont, E. Kestenare) 
CC                         add for POC in sediments
CC----------------------------------------------------------------------
CC parameters and commons
CC ======================

c#include "parameter.h"
c#include "common.h"
       USE myalloc

CC+CC  Implicit typing is never allowed
        IMPLICIT NONE
CC+CC  Implicit typing is never allowed

#if defined key_passivetrc
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER ji,jj,jk,jn
CC----------------------------------------------------------------------
CC statement functions
CC ===================

#include "stafun.h"

CCC---------------------------------------------------------------------
CCC  OPA8, LODYC (15/11/96)
CCC---------------------------------------------------------------------
C
C
C 0. initialisations
C ------------------
C
      IF(lwp) THEN
          WRITE(numout,*) ' '
          WRITE(numout,*) ' *** trcdtr initialisation for '
          WRITE(numout,*) '     passive tracers'
          WRITE(numout,*) ' '
      ENDIF

#if defined key_trc_generic
#include "trcdtr.generic.h"
#elif defined key_trc_npzd
#include "trcdtr.npzd.h"
#elif defined key_trc_npzdb
da inserire       
#elif key_trc_bfm
# include "trcdtr.bfm.h"

#endif
C
C before field :
C
       DO jn=1,jptra 
         DO  jk=1,jpk
           DO  jj=1,jpj
             DO  ji=1,jpi
               trb(ji,jj,jk,jn)=trn(ji,jj,jk,jn)
             END DO
           END DO
         END DO
       END DO

C

#else
C
C     no passive tracers
C
#endif
C
      RETURN
      END

