C $Id: trcsms.F,v 1.2 2009-09-11 09:20:56 cvsogs01 Exp $

      SUBROUTINE trcsms(ktask,kt)
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE trcsms
CCC                     *****************
CCC
CCC  PURPOSE :
CCC  ---------
CCC	time loop of opa for passive tracer
CCC
CC   METHOD :
CC   -------
CC      compute the well/spring evolution
CC
CC   INPUT :
CC   -----
CC      argument
CC              ktask           : task identificator
CC              kt              : time step
CC      COMMON
CC		all the COMMON defined in opa
CC
CC
CC   OUTPUT :			: no
CC   ------
CC
CC   WORKSPACE :
CC   ---------
CC
CC   EXTERNAL :
CC   --------
CC      trcbio, trcsed, trcopt for NPZD or NNPZDDOM models
CC
CC      h3cprg for HAMOC3 and P3ZD models
CC
CC
CC   MODIFICATIONS:
CC   --------------
CC      original  : 96-11
CC      additions : 99-07 (M. Levy)
CC                  04-00 (O. Aumont, M.A. Foujols) HAMOCC3 and P3ZD
CC                  12-00 (O. Aumont, E. Kestenare) add trcexp for 
CC                                                 instantaneous export
CC----------------------------------------------------------------------
CC parameters and commons
CC ======================

c#include "parameter.h"
c#include "common.h"
       USE myalloc

CC+CC  Implicit typing is never allowed
        IMPLICIT NONE
CC+CC  Implicit typing is never allowed
#include <mpif.h>

CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER ktask,kt
      INTEGER ji,jj,jk,jn

#if defined key_passivetrc
C
CC----------------------------------------------------------------------
CC statement functions
CC ===================

#include "stafun.h"

CCC---------------------------------------------------------------------
CCC  OPA8, LODYC (15/11/96)
CCC---------------------------------------------------------------------
C
C
C this ROUTINE is called only every ndttrc time step
C
CCC 10 11 2004  F79 cronometer-start

       trcsmsparttime = MPI_WTIME()

CCC

      IF (mod(kt,ndttrc) .NE. 0) RETURN
C
C this first routines are parallelized on vertical slab
C
C
#    if defined key_trc_npzd || defined key_trc_nnpzddom || key_trc_bfm
C
C tracers: optical model
C ----------------------
C
       CALL trcopt(ktask,kt)
C
C tracers: biological model
C -------------------------
C
       CALL trcbio(ktask,kt)
C
C tracers: sedimentation model
C ----------------------------
C
CCC trcsed-trexp no updated for time step advancing 
#if  defined key_trc_sed
       CALL trcsed(ktask,kt)
C       CALL trcexp
# endif
C
#    elif defined key_trc_hamocc3
C
C    h3cprg: main PROGRAM for HAMOCC3 or P3ZD models
C    -----------------------------------------------
       CALL h3cprg(ktask,kt)
#    else
C
C SMS to DO
C
#    endif
#else
C
C no passive tracers
C
#endif
C
CCC 10 11 2004  F79 cronometer-stop

       trcsmsparttime = MPI_WTIME() - trcsmsparttime
       trcsmstottime = trcsmstottime + trcsmsparttime

CC-CC       write(*,*) "F79T:trcsmsparttime", trcsmsparttime
CC-CC       write(*,*) "F79T:trcsmstottime", trcsmstottime

CCC


      RETURN
      END
