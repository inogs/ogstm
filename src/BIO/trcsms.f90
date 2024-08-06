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
       use simple_timer

       IMPLICIT NONE


!! this ROUTINE is called only every ndttrc time step

       trcsmsparttime = MPI_WTIME() ! cronometer-start

!! this first routines are parallelized on vertical slab

       call tstart("trcopt")

       CALL trcopt ! tracers: optical model

       call tstop("trcopt")

       call tstart("trcbio")
       CALL trcbio ! tracers: biological model
       call tstop("trcbio")

!! trcsed no updated for time step advancing
#if  defined key_trc_sed
       call tstart("trcsed")

       CALL trcsed ! tracers: sedimentation model

       call tstop("trcsed")
# endif

       trcsmsparttime = MPI_WTIME() - trcsmsparttime ! cronometer-stop
       trcsmstottime = trcsmstottime + trcsmsparttime
       
      END SUBROUTINE trcsms
