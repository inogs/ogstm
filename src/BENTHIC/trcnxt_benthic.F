      SUBROUTINE trcnxt_benthic

       USE myalloc
       USE myalloc_mpp
       USE BIOBEN_mem
       USE Time_Manager

      IMPLICIT NONE

CC local declarations
CC ==================
      INTEGER pack_size
! omp variables
      INTEGER :: mytid, ntids


#ifdef __OPENMP
      INTEGER ::  omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
      EXTERNAL :: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
#endif

      INTEGER ji,jj,jk,jn,jp


#ifdef __OPENMP
      ntids = omp_get_max_threads() ! take the number of threads
      mytid = -1000000
#else
      ntids = mpi_pack_size
      mytid = 0
#endif


       trcnxtparttime = MPI_WTIME() ! cronometer-start

C 1. fields at the next time
C --------------------------
C Tracer slab
C ===========


      TRACER_LOOP: DO  jn = 1, jptra_b, ntids


C 1. Lateral boundary conditions on tra (1,1,1,jn)

#ifdef key_mpp

C   ... Mpp : export boundary values to neighboring processors

        IF( ntids - 1 + jn <= jptra ) THEN
           pack_size = ntids
        ELSE
           pack_size = ntids - (ntids - 1 + jn - jptra)
        END IF

        CALL mpplnk_my(tra_ben(1,1,1,jn), pack_size,1,1)

#  else

C   ... T-point, 3D array, full array tra(1,1,1,jn) is initialised

        CALL lbc( tra_ben(1,1,1,jn), 1, 1, 1, 1, jpk, 1 )

#endif


!$omp   parallel default(none) private(mytid,ji,jj,jk)
!$omp&      shared(jn,jpk,jpj,jpi,trn,trb,tra,tmask) 

#ifdef __OPENMP
        mytid = omp_get_thread_num()  ! take the thread ID
#else
      PACK_LOOP1: DO jp=1,ntids
       mytid=jp-1
#endif
      IF( mytid + jn <= jptra_b ) THEN


        DO jk = 1,jpk_ben
          DO jj = 1,jpj
            DO ji = 1,jpi

            trn_ben(ji,jj,jk,jn+mytid) = trn_ben(ji,jj,jk,jn+mytid) + tra_ben(ji,jj,jk,jn+mytid)*rdt
            trb_ben(ji,jj,jk,jn+mytid) = trn_ben(ji,jj,jk,jn+mytid)
            tra_ben(ji,jj,jk,jn+mytid) = 0.e0

            END DO
          END DO
        END DO

      END IF
!$omp end parallel
#ifdef __OPENMP
#else
      END DO PACK_LOOP1
      mytid =0
#endif


C END of tracer slab
C ==================

       END DO TRACER_LOOP


       trcnxtparttime = MPI_WTIME() - trcnxtparttime ! cronometer-stop
       trcnxttottime = trcnxttottime + trcnxtparttime


      END SUBROUTINE trcnxt_benthic
