C $Id: h3cprg.F,v 1.2 2009-09-11 09:20:56 cvsogs01 Exp $

       SUBROUTINE h3cprg(ktask,kt)
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE h3cprg
CCC                     *****************
CCC
CCC  PURPOSE :
CCC  ---------
CCC     Call Biological sources and sinks subroutines
CCC
CC   METHOD :
CC   -------
CC      
CC
CC   INPUT :
CC   -----
CC      argument
CC              ktask           : task identificator
CC              kt              : time step
CC      common
CC              all the common defined in opa
CC
CC
CC   OUTPUT :                   : no
CC   ------
CC
CC   WORKSPACE :
CC   ---------
CC
CC   EXTERNAL :
CC   --------
CC      for HAMOCC3:
CC      h3cche, h3cint, h3clys, h3cbio, h3crem, h3cflx
CC      for O. Aumont P3ZD:
CC      h3cche, h3cint, p3zbio, h3cexp, p3zrem, h3cflx, 
CC
CC   MODIFICATIONS:
CC   --------------
CC      original  : E. Maier-Reimer (GBC 1993)
CC      additions : O. Aumont (1998)
CC      additions : C. Le Quere (1999)
CC      additions : O. Aumont (2001)
CC      additions : O. Aumont , EK (05/2001): add h3cadj 
CC----------------------------------------------------------------------
CC parameters and commons
CC ======================

c#include "parameter.h"
c#include "common.h"
       USE myalloc

CC+CC  Implicit typing is never allowed
        IMPLICIT NONE
CC+CC  Implicit typing is never allowed

CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER ktask, kt
#if defined key_passivetrc && defined key_trc_hamocc3
      INTEGER ji,jj,jk,jn
C
C this part is without macrotasking coding
C
      IF (ktask.EQ.1) THEN
C
C Compute tracers adjustments
C ---------------------------
C
          CALL h3cadj
C
C Compute chemical variables
C --------------------------
C
          CALL h3cche
C......................................................................
C
C Interpolate chemical variables
C ------------------------------
C
         CALL h3cint(kt)

C
C......................................................................
C
C Compute CaCO3 saturation
C ------------------------
C
          CALL h3clys
C
C......................................................................
C
C Compute soft tissue production (POC)
C ------------------------------------
C
#if defined key_trc_p3zd
          CALL p3zbio
#else
          CALL h3cbio
#endif
          CALL h3cexp
          CALL trcexp

C
C......................................................................
C
C Compute soft tissue remineralisation
C ------------------------------------
C
#if defined key_trc_p3zd
          CALL p3zrem
          trb(:,:,:,:)=trn(:,:,:,:)
#else
          CALL h3crem
#endif
C
C......................................................................
C
C Compute surface fluxes
C ----------------------
C
          CALL h3cflx
          CALL h3cvfl(kt)

      ENDIF 
C
C.............................synchro..................................
C
      CALL forkjoin
C
C......................................................................
C
#endif
C
      RETURN
      END

