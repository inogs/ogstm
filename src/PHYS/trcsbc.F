      SUBROUTINE trcsbc

! Computes surface boundary conditions on passive tracers

      USE myalloc
      USE myalloc_mpp

      IMPLICIT NONE

      INTEGER :: mytid, ntids


#ifdef __OPENMP
      INTEGER ::  omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
      EXTERNAL :: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
#endif


      INTEGER  :: ji,jj,jn,jb, diavar
      REAL(8)  :: ztra,zse3t

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

               zse3t = 1. / e3t(ji,jj,1)

                  ztra = 1./ rhopn(ji,jj,1) * zse3t * tmask(ji,jj,1) * emp(ji,jj) * trn(ji,jj,1,jn) ! original emps(ji,jj)
                  tra(ji,jj,1,jn) = tra(ji,jj,1,jn) + ztra

          END DO
        END DO
      ENDDO


      Do jn=ppO3c,ppO3h

         if ( jn==ppO3c )  diavar=ppO3c_concdilu_flux
         if ( jn==ppO3h )  diavar=ppO3h_concdilu_flux
         tra_DIA_2d(:,:,diavar) = 0.0
         DO jb=1,NBFMPOINTS_SUP,ntids


!$omp  parallel default(none) private(mytid,ji,jj,zse3t,ztra)
!$omp&      shared(jb,jn, rhopn,tmask,emp,trn,tra,tra_DIA_2d,NBFMPOINTS_SUP,BFMpoints,diavar)
#ifdef __OPENMP
        mytid = omp_get_thread_num()  ! take the thread ID
#endif
         IF( mytid + jb <= NBFMPOINTS_SUP) THEN
          ji = BFMpoints(1,jb+mytid)
          jj = BFMpoints(2,jb+mytid)

          ztra = 1./ rhopn(ji,jj,1) * zse3t * tmask(ji,jj,1) *emp(ji,jj) * trn(ji,jj,1,jn) ! original emps(ji,jj)

          tra_DIA_2d(ji,jj,diavar) = ztra ! save evaporation rate

         ENDIF
!$omp end parallel
         END DO
      ENDDO



         trcsbcparttime = MPI_WTIME()   - trcsbcparttime
         trcsbctottime  = trcsbctottime + trcsbcparttime
      END SUBROUTINE trcsbc
