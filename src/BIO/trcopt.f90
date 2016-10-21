
      SUBROUTINE trcopt
!!!---------------------------------------------------------------------
!!!
!!!                       ROUTINE trcopt
!!!                     ******************

       USE myalloc
       USE myalloc_mpp
       USE OPT_mem
       IMPLICIT NONE


!!! local declarations
!!! ==================

      REAL(8) conversion
#if defined key_trc_nnpzddom || defined key_trc_npzd || key_trc_bfm

      INTEGER :: jk,jj,ji
      INTEGER :: mytid, ntids! omp variables

#ifdef __OPENMP1
      INTEGER ::  omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
      EXTERNAL :: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
#endif


#ifdef __OPENMP1
      ntids = omp_get_max_threads() ! take the number of threads
      mytid = -1000000

#else
      ntids = 1
      mytid = 0
#endif



      trcoptparttime = MPI_WTIME() ! cronometer-start

      conversion = 0.50/0.217 ! Conversion Einstein to Watt  E2W=0.217

! vertical slab
! ===============

      DO jj = 1,jpj,ntids
!!!$omp parallel default(none) private(mytid,jk,ji)
!!!$omp&                       shared(jj,jpk,jpj,jpi,xpar,conversion,kef,e3t,qsr)
#ifdef __OPENMP1
         mytid  = omp_get_thread_num()  ! take the thread ID
#endif
      if (mytid+jj.le.jpj) then


! 1. determination of surface irradiance

        DO ji = 1,jpi

          zpar0m(ji)          = qsr(ji,mytid+jj)*conversion
          zpar100(ji)         = zpar0m(ji)*0.01
          xpar(ji,mytid+jj,1) = zpar0m(ji)
          zpar(ji,1)          = zpar0m(1)
          xEPS(ji,1)          = kef(ji,mytid+jj)

        END DO

!! 2. determination of xpar
!! ------------------------

        DO jk = 2,jpk
          DO ji = 1,jpi

            xEPS(jk,ji)          = kef(ji,mytid+jj)
            xEPS(jk,ji)          = max(xEPS(jk,ji),1.D-15) ! avoid denormalized number
            xpar(jk,mytid+jj,ji) = xpar(jk,mytid+jj,ji-1) *exp(-1. * xEPS(jk-1,ji)* e3t(jk-1,jj,ji))
            xpar(jk,mytid+jj,ji) = max(xpar(jk,mytid+jj,ji),1.D-15) ! avoid denormalized number

          END DO
        END DO

        DO jk = 1,jpk
          DO ji = 1,jpi
            !a=xpar(jk,mytid+jj,ji)  xpar(jk,mytid+jj,ji) = max(a*exp(- xEPS(jk,ji)* 0.5D+00* e3t(jk) ), 1.D-15)
      
            xpar(jk,mytid+jj,ji) = xpar(jk,mytid+jj,ji) * exp(- xEPS(jk,ji)* 0.5D+00* e3t(jk,jj,ji) )
            xpar(jk,mytid+jj,ji) = max(xpar(jk,mytid+jj,ji),1.D-15)
          END DO
        END DO

      endif
!!!$omp end parallel
      END DO

       trcoptparttime = MPI_WTIME() - trcoptparttime ! cronometer-stop
       trcopttottime = trcopttottime + trcoptparttime

#else

!!    No optical model

#endif

      END SUBROUTINE trcopt
