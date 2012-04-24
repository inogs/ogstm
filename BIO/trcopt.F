
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

#if defined key_trc_nnpzddom || defined key_trc_npzd || key_trc_bfm

      INTEGER ji,jj,jk,jn

       trcoptparttime = MPI_WTIME() ! cronometer-start

C vertical slab
C ===============

      DO jj = 1,jpj

C 1. determination of surface irradiance
C
        DO ji = 1,jpi

          zpar0m(ji)    = qsr(ji,jj)*0.50/0.217 ! Conversion Einstein to Watt  E2W=0.217
          zpar100(ji)   = zpar0m(ji)*0.01
          xpar(ji,jj,1) = zpar0m(ji)
          zpar(ji,1)    = zpar0m(1)

        END DO
C
C 2. determination of xpar
C ------------------------
C
        DO jk = 2,jpk
          DO ji = 1,jpi

            xEPS(ji,jk) = kef(ji,jj)

            xEPS(ji,jk) = max(xEPS(ji,jk)
     &          ,1.D-15) ! avoid denormalized number

             xpar(ji,jj,jk) = xpar(ji,jj,jk-1) *
     &                    exp(-1. * xEPS(ji,jk-1)* e3t(jk-1))


            xpar(ji,jj,jk) = max(xpar(ji,jj,jk)
     &          ,1.D-15) ! avoid denormalized number

          END DO
        END DO

        DO jk = 1,jpk
          DO ji = 1,jpi
            xpar(ji,jj,jk) = xpar(ji,jj,jk) * exp(- xEPS(ji,jk)* 
     &      0.5D+00* e3t(jk) )

            xpar(ji,jj,jk) = max(xpar(ji,jj,jk)
     &          ,1.D-15)
          END DO
        END DO

      END DO

       trcoptparttime = MPI_WTIME() - trcoptparttime !F79 cronometer-stop
       trcopttottime = trcopttottime + trcoptparttime

#    else

C    No optical model

#    endif

      END SUBROUTINE trcopt
