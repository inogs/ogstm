      SUBROUTINE trcsms
!!!---------------------------------------------------------------------
!!!
!!!                       ROUTINE trcsms
!!!                     *****************
!!!
!!!  PURPOSE :
!!!  ---------
!!!    time loop of ogstm for passive tracer
!!!
!!   METHOD :
!!   -------
!!      compute the well/spring evolution
!!
!!   INPUT :
!!   -----
!!
!!   EXTERNAL :
!!   --------
!!      trcbio, trcsed, trcopt for NPZD or NNPZDDOM models


       USE myalloc
       USE mpi
       IMPLICIT NONE


!! this ROUTINE is called only every ndttrc time step

       trcsmsparttime = MPI_WTIME() ! cronometer-start


!! this first routines are parallelized on vertical slab

       CALL trcopt ! tracers: optical model

       CALL trc3streams
       
       CALL trcbio ! tracers: biological model

!! trcsed no updated for time step advancing
#if  defined key_trc_sed
       CALL trcsed ! tracers: sedimentation model
# endif

       trcsmsparttime = MPI_WTIME() - trcsmsparttime ! cronometer-stop
       trcsmstottime = trcsmstottime + trcsmsparttime
       
      END SUBROUTINE trcsms
