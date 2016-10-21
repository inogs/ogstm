      SUBROUTINE trcsbc

! Computes surface boundary conditions on passive tracers

      USE myalloc
      USE myalloc_mpp

      IMPLICIT NONE

      INTEGER :: mytid, ntids


#ifdef __OPENMP1
      INTEGER ::  omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
      EXTERNAL :: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
#endif


      INTEGER  :: jj,jijn
      REAL(8)  :: ztra,zse3t

#ifdef __OPENMP1
      ntids = omp_get_max_threads() ! take the number of threads
      mytid = -1000000
#else
      ntids = 1
      mytid = 0
#endif


      trcsbcparttime = MPI_WTIME()

! Conc/dilution process


      Do jn=1,jptra
        DO jj = 1, jpj
            DO ji = 1, jpi

               zse3t = 1. / e3t(1,jj,ji)

                  ztra = 1./ rhopn(1,jj,ji) * zse3t * tmask(1,jj,ji) * emp(jj,ji) * trn(1,jj,ji,jn) ! original emps(jj,ji)
                  tra(1,jj,ji,jn) = tra(1,jj,ji,jn) + ztra

          END DO
        END DO
      ENDDO

         trcsbcparttime = MPI_WTIME()   - trcsbcparttime
         trcsbctottime  = trcsbctottime + trcsbcparttime
      END SUBROUTINE trcsbc
