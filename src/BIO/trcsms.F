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
       USE myalloc_mpp
       IMPLICIT NONE


C this ROUTINE is called only every ndttrc time step

       trcsmsparttime = MPI_WTIME() ! cronometer-start


C this first routines are parallelized on vertical slab

!      CALL trcopt ! tracers: optical model
       CALL trcbio ! tracers: biological model

C trcsed no updated for time step advancing
#if  defined key_trc_sed
       CALL trcsed ! tracers: sedimentation model
# endif

       trcsmsparttime = MPI_WTIME() - trcsmsparttime ! cronometer-stop
       trcsmstottime = trcsmstottime + trcsmsparttime

      END SUBROUTINE trcsms
