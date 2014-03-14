      SUBROUTINE trcsbc

! Computes surface boundary conditions on passive tracers

      USE myalloc
      USE stafun

      IMPLICIT NONE

      INTEGER :: mytid, ntids, itid


#ifdef __OPENMP
      INTEGER ::  omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
      EXTERNAL :: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
#endif


      INTEGER  :: ji,jj,jn,jb
      REAL(8)  :: ztra,zse3t
#include <mpif.h>
!      zsrau = 1. / rau0

#ifdef __OPENMP
      ntids = omp_get_max_threads() ! take the number of threads
      mytid = -1000000
#else
      ntids = 1
      mytid = 0
#endif


      trcsbcparttime = MPI_WTIME()

! Conc/dilution process

      Do jn=50,jptra
        DO jj = 1, jpj
            DO ji = 1, jpi

                  zse3t = 1. / fse3t(ji,jj,1)

                  ztra = 1./ rhopn(ji,jj,1) * zse3t * tmask(ji,jj,1) * emp(ji,jj) * trn(ji,jj,1,jn) ! original emps(ji,jj)
                  tra(ji,jj,1,jn) = tra(ji,jj,1,jn) + ztra
                  tra_pp(ji,jj,jn-49,jptra_dia) = ztra ! save evaporation rate

          END DO
        END DO
      ENDDO



      Do jn=50,jptra
      tra_pp(:,:,jn-49,jptra_dia) = 0.0

         DO jb=1,NBFMPOINTS_SUP,ntids

!$omp  parallel default(none) private(mytid,ji,jj,zse3t,ztra)
!$omp&      shared(jb,jn, rhopn,tmask,emp,trn,tra, tra_pp,NBFMPOINTS_SUP,BFMpoints)
#ifdef __OPENMP
        mytid = omp_get_thread_num()  ! take the thread ID
#endif
         IF( mytid + jb <= NBFMPOINTS_SUP) THEN
          ji = BFMpoints(1,jb+mytid)
          jj = BFMpoints(2,jb+mytid)

          zse3t = 1. / fse3t(ji,jj,1)

          ztra = 1./ rhopn(ji,jj,1) * zse3t * tmask(ji,jj,1) * emp(ji,jj) * trn(ji,jj,1,jn) ! original emps(ji,jj)

          tra_pp(ji,jj,jn-49,jptra_dia) = ztra ! save evaporation rate

         ENDIF
!$omp end parallel
         END DO
      ENDDO






! Contribution of surface elevation virtual flux

         trcsbcparttime = MPI_WTIME()   - trcsbcparttime
         trcsbctottime  = trcsbctottime + trcsbcparttime
      END SUBROUTINE trcsbc
