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

       ! XXX: to remove
       use BIO_mem, only: ogstm_sediPI,ogstm_PH,ogstm_co2
       USE OPT_mem, only: kef
       IMPLICIT NONE


!! this ROUTINE is called only every ndttrc time step

       trcsmsparttime = MPI_WTIME() ! cronometer-start


!! this first routines are parallelized on vertical slab

       call tstart("trcopt")

       !$acc update device(kef,qsr)
       CALL trcopt ! tracers: optical model

       call tstop("trcopt")

       call tstart("trcbio")
       !$acc update device(mbathy,bfmmask,trn,DAY_LENGTH,vatm,tn,sn,rho,e3t,gdept,ogstm_PH,ogstm_co2)
       CALL trcbio ! tracers: biological model
       !$acc update host(tra,tra_DIA,tra_DIA_2d,ogstm_sediPI,ogstm_PH)
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
