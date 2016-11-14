      SUBROUTINE trcsbc

! Computes surface boundary conditions on passive tracers

      USE myalloc
      USE mpi
      ! epascolo USE myalloc_mpp

      IMPLICIT NONE

      INTEGER :: mytid, ntids

      INTEGER  :: jj,ji,jn
      double precision  :: ztra,zse3t

      trcsbcparttime = MPI_WTIME()

! Conc/dilution process


                  DO jn=1,jptra
            DO ji = 1, jpi
        DO jj = 1, jpj

               zse3t = 1. / e3t(1,jj,ji)

                  ztra = 1./ rhopn(1,jj,ji) * zse3t * tmask(1,jj,ji) * emp(jj,ji) * trn(1,jj,ji,jn) ! original emps(jj,ji)
                  tra(1,jj,ji,jn) = tra(1,jj,ji,jn) + ztra

          END DO
        END DO
      ENDDO

         trcsbcparttime = MPI_WTIME()   - trcsbcparttime
         trcsbctottime  = trcsbctottime + trcsbcparttime
      END SUBROUTINE trcsbc
