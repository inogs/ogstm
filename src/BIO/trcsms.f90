      SUBROUTINE trcsms(datestring)
!!!---------------------------------------------------------------------
!!!
!!!                       ROUTINE trcsms
!!!                     *****************
!!!
!!!  PURPOSE :
!!!  ---------
!!!
!!   METHOD :
!!   -------
!!
!!   INPUT :
!!   -----
!!
!!   EXTERNAL :
!!   --------


       USE myalloc
       USE mpi
       IMPLICIT NONE

       character(LEN=17), INTENT(IN) ::  datestring


!! this ROUTINE is called only every ndttrc time step

       trcsmsparttime = MPI_WTIME() ! cronometer-start


!! this first routines are parallelized on vertical slab

       CALL trcopt ! tracers: optical model

       CALL trc3streams(datestring)
       
       CALL trcbio ! tracers: biological model

!! trcsed no updated for time step advancing
#if  defined key_trc_sed
       CALL trcsed ! tracers: sedimentation model
# endif

       trcsmsparttime = MPI_WTIME() - trcsmsparttime ! cronometer-stop
       trcsmstottime = trcsmstottime + trcsmsparttime
       
      END SUBROUTINE trcsms
