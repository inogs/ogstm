      SUBROUTINE trcsbc

! Computes surface boundary conditions on passive tracers

      USE myalloc
      USE stafun

      IMPLICIT NONE

      INTEGER :: mytid, ntids


#ifdef __OPENMP
      INTEGER ::  omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
      EXTERNAL :: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
#endif


      INTEGER  :: ji,jj,jn
      REAL(8)  :: ztra,zse3t
#include <mpif.h>

#ifdef __OPENMP
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

                  zse3t = 1. / fse3t(ji,jj,1)

                  ztra = 1./ rhopn(ji,jj,1) * zse3t * tmask(ji,jj,1) * emp(ji,jj) * trn(ji,jj,1,jn) ! original emps(ji,jj)
                  tra(ji,jj,1,jn) = tra(ji,jj,1,jn) + ztra

          END DO
        END DO
      ENDDO

         trcsbcparttime = MPI_WTIME()   - trcsbcparttime
         trcsbctottime  = trcsbctottime + trcsbcparttime
      END SUBROUTINE trcsbc
