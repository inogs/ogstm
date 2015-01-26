
      SUBROUTINE trcopt(ktask,kt)
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE trcopt
CCC                     ******************

       USE myalloc
       USE myalloc_mpp
       USE OPT_mem
       IMPLICIT NONE


CC local declarations
CC ==================
      INTEGER kt,ktask
      REAL(8) conversion
#if defined key_trc_nnpzddom || defined key_trc_npzd || key_trc_bfm

      INTEGER :: ji,jj,jk,jn
      INTEGER :: mytid, ntids! omp variables

#ifdef __OPENMP
      INTEGER ::  omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
      EXTERNAL :: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
#endif


#ifdef __OPENMP
      ntids = omp_get_max_threads() ! take the number of threads
      mytid = -1000000

#else
      ntids = 1
      mytid = 0
#endif



      trcoptparttime = MPI_WTIME() ! cronometer-start

      conversion = 0.50/0.217 ! Conversion Einstein to Watt  E2W=0.217

C vertical slab
C ===============

      DO jj = 1,jpj,ntids
!$omp parallel default(none) private(mytid,ji,jk)
!$omp&                       shared(jj,jpi,jpj,jpk,xpar,conversion,kef,e3t,qsr)
#ifdef __OPENMP
         mytid  = omp_get_thread_num()  ! take the thread ID
#endif
      if (mytid+jj.le.jpj) then


C 1. determination of surface irradiance
C
        DO ji = 1,jpi

          zpar0m(ji)          = qsr(ji,mytid+jj)*conversion
          zpar100(ji)         = zpar0m(ji)*0.01
          xpar(ji,mytid+jj,1) = zpar0m(ji)
          zpar(ji,1)          = zpar0m(1)
          xEPS(ji,1)          = kef(ji,mytid+jj)

        END DO
C
C 2. determination of xpar
C ------------------------
C
        DO jk = 2,jpk
          DO ji = 1,jpi

            xEPS(ji,jk)          = kef(ji,mytid+jj)
            xEPS(ji,jk)          = max(xEPS(ji,jk),1.D-15) ! avoid denormalized number
            xpar(ji,mytid+jj,jk) = xpar(ji,mytid+jj,jk-1) *exp(-1. * xEPS(ji,jk-1)* e3t(jk-1))
            xpar(ji,mytid+jj,jk) = max(xpar(ji,mytid+jj,jk),1.D-15) ! avoid denormalized number

          END DO
        END DO

        DO jk = 1,jpk
          DO ji = 1,jpi
            !a=xpar(ji,mytid+jj,jk); xpar(ji,mytid+jj,jk) = max(a*exp(- xEPS(ji,jk)* 0.5D+00* e3t(jk) ), 1.D-15);
            xpar(ji,mytid+jj,jk) = xpar(ji,mytid+jj,jk) * exp(- xEPS(ji,jk)* 0.5D+00* e3t(jk) )
            xpar(ji,mytid+jj,jk) = max(xpar(ji,mytid+jj,jk),1.D-15)
          END DO
        END DO

      endif
!$omp end parallel
      END DO

       trcoptparttime = MPI_WTIME() - trcoptparttime !F79 cronometer-stop
       trcopttottime = trcopttottime + trcoptparttime

#    else

C    No optical model

#    endif

      END SUBROUTINE trcopt
