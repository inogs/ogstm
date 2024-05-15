
      SUBROUTINE trcopt
!!!---------------------------------------------------------------------
!!!
!!!                       ROUTINE trcopt
!!!                     ******************

       USE myalloc
       USE mpi
       ! epascolo USE myalloc_mpp
       USE OPT_mem
       IMPLICIT NONE


!!! local declarations
!!! ==================

      double precision :: conversion
#if defined key_trc_nnpzddom || defined key_trc_npzd || key_trc_bfm

      INTEGER :: jk,jj,ji,queue
      INTEGER :: mytid, ntids! omp variables

      trcoptparttime = MPI_WTIME() ! cronometer-start

      conversion = 0.50/0.217 ! Conversion Einstein to Watt  E2W=0.217

! vertical slab
! ===============

      
!!!$omp parallel default(none) private(mytid,jk,ji)
!!!$omp&                       shared(jj,jpk,jpj,jpi,xpar,conversion,kef,e3t,qsr)

! 1. determination of surface irradiance
      queue=1
      ! XXX: this can be improved but no way to test as it doesn't seem to
      ! impact numerical results
      !$acc kernels default(present) async(queue)
      DO ji = 1,jpi
         DO jj = 1,jpj

            ! zpar0m(ji)          = qsr(jj,ji)*conversion
            ! zpar100(ji)         = zpar0m(ji)*0.01
            xpar(1,jj,ji) = qsr(jj,ji)*conversion
            ! zpar(1,ji)          = zpar0m(1)
            xEPS_ogstm(1,ji)          = kef(jj,ji)

            !! 2. determination of xpar
            !! ------------------------

            DO jk = 2,jpk
               xEPS_ogstm(jk,ji) = max(kef(jj,ji),1.D-15) ! avoid denormalized number
            END DO

            DO jk = 2,jpk
               xpar(jk,jj,ji) = max( xpar(jk-1,jj,ji) *exp(-1. * xEPS_ogstm(jk-1,ji)* e3t(jk-1,jj,ji)) ,1.D-15) ! avoid denormalized number

            END DO

            DO jk = 1,jpk
               xpar(jk,jj,ji) = max( xpar(jk,jj,ji) * exp(- xEPS_ogstm(jk,ji)* 0.5D+00* e3t(jk,jj,ji) ) ,1.D-15)
            END DO

         END DO
      ENDDO
      !$acc end kernels
      !$acc wait(queue)

       trcoptparttime = MPI_WTIME() - trcoptparttime ! cronometer-stop
       trcopttottime = trcopttottime + trcoptparttime

#else

!!    No optical model

#endif

      END SUBROUTINE trcopt
