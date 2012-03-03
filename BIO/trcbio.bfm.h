
       USE myalloc
       USE myalloc_mpp
       USE BIO_mem
       USE BC_mem

       IMPLICIT NONE


CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER kt,ktask
      LOGICAL sur,bot,isT,isBIO
      REAL(8) a(jptra),b(jptra),c(4),d(jptra_dia)

      INTEGER ji,jj,jk,jn
      INTEGER jtr,jtrmax,tra_idx
      REAL(8):: opa_den, opa_ice, opa_co2
! omp variables
            INTEGER :: mytid, ntids, itid

#ifdef __OPENMP
            INTEGER ::  omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
            EXTERNAL :: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
#endif

CC----------------------------------------------------------------------
CC statement functions
CC ===================

#include "stafun.h"

      INTERFACE OPA_Output_EcologyDynamics
         subroutine OPA_Output_EcologyDynamics(opa_tra, dim_opa_tra, sediPI_P1, 
     &       sediPI_P2, sediPI_P3,  sediPI_P4, local_opa_dia)
!            use global_mem, ONLY:RLEN
            IMPLICIT NONE
            integer dim_opa_tra
            real(8):: sediPI_P1, sediPI_P2, sediPI_P3, sediPI_P4
            real(8):: opa_tra(dim_opa_tra)
            real(8):: local_opa_dia(10)
         end subroutine
      END INTERFACE

C   | --------------|
C   | BFM MODEL CALL|
C   | --------------|
C
#ifdef __OPENMP
      ntids = omp_get_max_threads() ! take the number of threads
      mytid = -1000000
#else
      ntids = 1
      mytid = 0
#endif
          surf_mask(:) = 0.
          surf_mask(1) = 1.

          opa_den=1029
          opa_ice=0
          opa_co2=365.0
          tra_idx = tra_matrix_gib(1)
          jtrmax=jptra

      MAIN_LOOP: DO  jj = 1, jpjm1, ntids

!$omp   parallel default(none) private(jk,ji,mytid,isT,isBIO,sur,bot,jtr,a,b,d)
!$omp&      shared(jj,jpjm1,jpkbm1,jpim1,Tmask,tra_idx,tra_matrix_gib,
!$omp&		       restotr,jtrmax,trn,tn,sn,xpar,e3t,vatm,surf_mask,
!$omp&             sediPI,tra_pp,tra,rhopn,opa_ice,opa_co2)

#ifdef __OPENMP
        mytid = omp_get_thread_num()  ! take the thread ID
#endif

	             IF( mytid + jj <= jpjm1 ) THEN
C
C 1. biological level
C ===================
C
                 DO jk=1,jpkbm1
                    DO ji = 2,jpim1

                       isT     = Tmask(ji,jj+mytid,jk) .eq. 1
                       isBIO   = restotr(ji,jj+mytid,jk,tra_idx) .eq. 0 ! no nudging points are considered

                       IF( isT .and. isBIO) THEN

                          IF(jk .eq. 1) sur = .TRUE.
                          bot = .FALSE.

                          DO jtr=1, jtrmax
                             a(jtr) = trn(ji,jj+mytid,jk,jtr)
                          END DO

                          call OPA_Input_EcologyDynamics(a,jtrmax,
     &                      tn(ji,jj+mytid,jk), sn(ji,jj+mytid,jk),
     &                      rhopn(ji,jj+mytid,jk), opa_ice, opa_co2, xpar(ji,jj+mytid,jk),
     &                      e3t(jk), sur, vatm(ji,jj+mytid) * surf_mask(jk),bot )

                          call OPA_reset()

                          call EcologyDynamics()

                          call OPA_Output_EcologyDynamics(b, jtrmax, sediPI(ji,jj+mytid,jk,1),
     &                          sediPI(ji,jj+mytid,jk,2), sediPI(ji,jj+mytid,jk,3),
     &                          sediPI(ji,jj+mytid,jk,4),d)

C                           write(*,*) 'sediPI -->ji',ji,'sediPI -->jj' ,jj+mytid
C                           write(*,*) 'sediPI -->jk',jk
C                           write(*,*) 'sediPI(ji,jj+mytid,jk,1) -->', sediPI(ji,jj+mytid,jk,1)
C                           write(*,*) 'sediPI(ji,jj+mytid,jk,2) -->', sediPI(ji,jj+mytid,jk,2)
C                           write(*,*) 'sediPI(ji,jj+mytid,jk,3) -->', sediPI(ji,jj+mytid,jk,3)
C                           write(*,*) 'sediPI(ji,jj+mytid,jk,4) -->', sediPI(ji,jj+mytid,jk,4)

                          DO jtr=1, jtrmax
                             tra(ji,jj+mytid,jk,jtr) =tra(ji,jj+mytid,jk,jtr) +b(jtr) ! trend
                          END DO

                          DO jtr=1,jptra_dia
                             tra_pp(ji,jj+mytid,jk,jtr) = d(jtr) ! diagnostic
                          END DO

                       ELSE
                          sediPI(ji,jj+mytid,jk,:)=0
                          tra_pp(ji,jj+mytid,jk,:)=0
                       ENDIF

                    END DO

                END DO

             ENDIF

!$omp end parallel

C END of slab
C ===========
C
                END DO MAIN_LOOP

