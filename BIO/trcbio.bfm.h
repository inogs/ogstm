c $Id: trcbio.bfm.h,v 1.2 2009-09-11 09:20:56 cvsogs01 Exp $
CCC
CCC      trcbio.npzd.h
CCC      *****************
CCC
CC   defined key : key_trc_npzd
CC   ===========
CC
CC
CC   INPUT :
CC   -----
CC      argument
CC              ktask           : task identificator
CC              kt              : time step
CC      COMMON
CC            /comcoo/          : orthogonal curvilinear coordinates
CC                                and scale factors
CC                                depths
CC            /cottrp/          : present and next fields for passive
CC                              : tracers
CC            /comtsk/          : multitasking
CC            /cotbio/          : biological parameters
CC
CC   OUTPUT :
CC   ------
CC      COMMON
CC            /cottrp/ tra      : general tracer trend increased by the
CC                                now horizontal tracer advection trend
CC            /cottbd/ trbio    : now horizontal tracer advection trend
CC                                (IF 'key_trc_diabio' is activated)
CC
CC   WORKSPACE :
CC   ---------
CC      local
CC               zdet,zzoo,zphy,znut              : now concentrations
CC               zlt,zlnut,zlpe                   : limitation terms for phyto
CC                                                 
CC               zflxnp,zflxpn,zflxpz,zflxdz      : fluxes between bio
CC                                                  boxes
CC               zflxpd,zflxzd,zflxdn
CC               zphya,zzooa,znuta,zdeta          : after bio trends
CC               zphimp, zmp, zphimz, zmz         : mortality terms
CC               zppz, zpdz, zpppz, zppdz, zfood  : preferences terms
CC               zfilpz, zfilpd                   : filtration terms
CC
CC   EXTERNAL :                   no
CC   --------
CC
CC   REFERENCES :                 no
CC   ----------
CC
CC   MODIFICATIONS:
CC   --------------
CC       original : 95-02 (M. Levy)
CC                  99-07 (M. Levy) version .h
CC                  99-09 (M. Levy) version with no deep mixing limitation
CC       adaptations : 00-12 (E. Kestenare) 
CC                     assign a parameter to name individual tracers
CC                     01-02 (E. Kestenare)
CC                     introduce jpno3 instead of jpnut 
CC                     01-02 (E. Kestenare) add sediments 
CC----------------------------------------------------------------------
       USE myalloc
       USE myalloc_mpp
       USE BIO_mem

CC+CC  Implicit typing is never allowed
        IMPLICIT NONE
CC+CC  Implicit typing is never allowed

CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER kt,ktask
      LOGICAL sur,bot
      REAL(8) a(jptra),b(jptra)

      INTEGER ji,jj,jk,jn
      INTEGER jtr,jtrmax
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
     &       sediPI_P2, sediPI_P3,  sediPI_P4, local_opa_ppg, local_opa_ppn, local_opa_ppb)
!            use global_mem, ONLY:RLEN
            IMPLICIT NONE
            integer dim_opa_tra
            real(8):: sediPI_P1, sediPI_P2, sediPI_P3, sediPI_P4
            real(8):: opa_tra(dim_opa_tra)
            real(8):: local_opa_ppg, local_opa_ppn, local_opa_ppb
         end subroutine
      END INTERFACE

ccC---------------------------------------------------------------------
CCC  OPA8, LODYC (15/11/96)
CCC---------------------------------------------------------------------
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
      jtrmax=jptra

      MAIN_LOOP: DO  jj = 1, jpjm1, ntids

!$omp   parallel default(none) private(jk,ji,mytid,sur,bot,jtr,a,b)
!$omp&      shared(jj,jpjm1,jpkbm1,jpim1,Tmask,jtrmax,trn,tn,sn,xpar,e3t,vatm,surf_mask,
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

                       IF(Tmask(ji,jj+mytid,jk) .eq. 1) THEN

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
     &                          sediPI(ji,jj+mytid,jk,4),tra_pp(ji,jj+mytid,jk,1),
     &                          tra_pp(ji,jj+mytid,jk,2),tra_pp(ji,jj+mytid,jk,3))

C                           write(*,*) 'sediPI -->ji',ji,'sediPI -->jj' ,jj+mytid
C                           write(*,*) 'sediPI -->jk',jk
C                           write(*,*) 'sediPI(ji,jj+mytid,jk,1) -->', sediPI(ji,jj+mytid,jk,1)
C                           write(*,*) 'sediPI(ji,jj+mytid,jk,2) -->', sediPI(ji,jj+mytid,jk,2)
C                           write(*,*) 'sediPI(ji,jj+mytid,jk,3) -->', sediPI(ji,jj+mytid,jk,3)
C                           write(*,*) 'sediPI(ji,jj+mytid,jk,4) -->', sediPI(ji,jj+mytid,jk,4)

                          DO jtr=1, jtrmax
                             tra(ji,jj+mytid,jk,jtr) =tra(ji,jj+mytid,jk,jtr) +b(jtr)
                          END DO

                       ENDIF

                    END DO

                END DO

             ENDIF

!$omp end parallel

C END of slab
C ===========
C
                END DO MAIN_LOOP

