
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
      !$acc parallel loop gang vector default(present) collapse(3) async(queue)
      DO ji = 1,jpi
         DO jj = 1,jpj
            DO jk = 1,jpk
               if (jk .eq. 1) then
                  xEPS_ogstm(jk,jj,ji)          = kef(jj,ji)
               else
                  xEPS_ogstm(jk,jj,ji) = max(kef(jj,ji),1.D-15) ! avoid denormalized number
               endif
            END DO
         END DO
      END DO

      !! 2. determination of xpar
      !! ------------------------

      !$acc parallel loop gang vector default(present) collapse(2) async(queue)
      ! NOTE: only collapse the first two loops as the third isn't parallel
      ! because we read xpar(jk-1,...) and write on xpar(jk,...)
      DO ji = 1,jpi
         DO jj = 1,jpj
            DO jk = 1,jpk
               if (jk .eq. 1) then
                  xpar(jk,jj,ji) = qsr(jj,ji)*conversion
               else
                  xpar(jk,jj,ji) = max( xpar(jk-1,jj,ji) *exp(-1. * xEPS_ogstm(jk-1,jj,ji)* e3t(jk-1,jj,ji)) ,1.D-15) ! avoid denormalized number
               endif
            END DO
         END DO
      END DO

      !$acc parallel loop gang vector default(present) collapse(3) async(queue)
      DO ji = 1,jpi
         DO jj = 1,jpj
            DO jk = 1,jpk
               xpar(jk,jj,ji) = max( xpar(jk,jj,ji) * exp(- xEPS_ogstm(jk,jj,ji)* 0.5D+00* e3t(jk,jj,ji) ) ,1.D-15)
            END DO
         END DO
      ENDDO
      !$acc wait(queue)

       trcoptparttime = MPI_WTIME() - trcoptparttime ! cronometer-stop
       trcopttottime = trcopttottime + trcoptparttime

#else

!!    No optical model

#endif

      END SUBROUTINE trcopt
