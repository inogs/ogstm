      SUBROUTINE trcsms(ktask,kt)
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE trcsms
CCC                     *****************
CCC
CCC  PURPOSE :
CCC  ---------
CCC    time loop of ogstm for passive tracer
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

CC
CC   EXTERNAL :
CC   --------
CC      trcbio, trcsed, trcopt for NPZD or NNPZDDOM models
CC
CC      h3cprg for HAMOC3 and P3ZD models

       USE myalloc
       IMPLICIT NONE

#include <mpif.h>
      INTEGER ktask,kt

#if defined key_passivetrc

C this ROUTINE is called only every ndttrc time step

       trcsmsparttime = MPI_WTIME() ! F79 cronometer-start

      IF (mod(kt,ndttrc) .NE. 0) RETURN

C this first routines are parallelized on vertical slab

#    if defined key_trc_npzd || defined key_trc_nnpzddom || key_trc_bfm

       CALL trcopt(ktask,kt) !tracers: optical model
       CALL trcbio(ktask,kt) ! tracers: biological model

CCC trcsed-trexp no updated for time step advancing 
#if  defined key_trc_sed
       CALL trcsed(ktask,kt) ! tracers: sedimentation model
C       CALL trcexp
# endif

#    elif defined key_trc_hamocc3

C     h3cprg: main PROGRAM for HAMOCC3 or P3ZD models
C    -----------------------------------------------
       CALL h3cprg(ktask,kt)
#    else

C SMS to DO

#    endif
#else

C no passive tracers

#endif


       trcsmsparttime = MPI_WTIME() - trcsmsparttime ! F79 cronometer-stop
       trcsmstottime = trcsmstottime + trcsmsparttime

      END SUBROUTINE trcsms
