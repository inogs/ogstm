      SUBROUTINE trcsms
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
CC
CC   EXTERNAL :
CC   --------
CC      trcbio, trcsed, trcopt for NPZD or NNPZDDOM models


       USE myalloc
       IMPLICIT NONE

#include <mpif.h>

#if defined key_passivetrc

C this ROUTINE is called only every ndttrc time step

       trcsmsparttime = MPI_WTIME() ! F79 cronometer-start


C this first routines are parallelized on vertical slab

#    if defined key_trc_npzd || defined key_trc_nnpzddom || key_trc_bfm

       CALL trcopt ! tracers: optical model
       CALL trcbio ! tracers: biological model

CCC trcsed-trexp no updated for time step advancing 
#if  defined key_trc_sed
       CALL trcsed ! tracers: sedimentation model
C      CALL trcexp
# endif

#    elif defined key_trc_hamocc3
#    else
#    endif
#else

C no passive tracers

#endif


       trcsmsparttime = MPI_WTIME() - trcsmsparttime ! F79 cronometer-stop
       trcsmstottime = trcsmstottime + trcsmsparttime

      END SUBROUTINE trcsms
