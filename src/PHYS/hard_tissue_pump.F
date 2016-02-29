
      SUBROUTINE hard_tissue_pump()
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE hard_tissue_pump
CCC                     *******************
CCC
CCC  PURPOSE :
CCC  ---------
CCC	compensate along the water column to preserve total mass
CCC
CCC
CC   METHOD :
CC   -------
CC    Must be called between trczdf and trcnxt
CC
C
CC
CC   OUTPUT :
CC   ------
CC
CC
CC   EXTERNAL :
CC   --------
CC
CC   REFERENCES :
CC   ----------
CC
CC   MODIFICATIONS:
CC   --------------

CC----------------------------------------------------------------------

       USE myalloc
       USE myalloc_mpp
       USE FN_mem
       USE BIO_mem

       IMPLICIT NONE


CC----------------------------------------------------------------------
CC local declarations
CC ==================

! omp variables
      INTEGER :: mytid, ntids, itid

#ifdef __OPENMP
      INTEGER ::  omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
      EXTERNAL :: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
#endif

      INTEGER ji,jj,jk,jn,jv
      INTEGER jpe,jpf
      REAL(8) RR_P2
      REAL(8) NORM_1, dCa
      REAL(8) g0,g1,dDic,dAlk
      REAL(8) SinkCa(jpk)


CC----------------------------------------------------------------------
CC statement functions
CC ===================


#ifdef __OPENMP
      ntids = omp_get_max_threads() ! take the number of threads
      mytid = -1000000
#else
      ntids = 1
      mytid = 0
#endif

! ***************************************************
      jpe    = 26
      jpf    = 2 ! CaCO3 dissolution along the whole water column  
      RR_P2  = 0.2   ![]
      dCa    = 1500. ![m]   Mediterranean Sea mean depth
      SinkCa = 0.
C Compute vertical integral
      DO jv  = 1,dimen_jvsnu, ntids ! cicla sui punti terra 2D

!$omp   parallel default(none) private(mytid,ji,jj,jk,jn,g1,SinkCa,dDic,dAlk)
!$omp&  shared(dimen_jvsnu,jarr_snu,jv,jpk,tmask,TOTcalc,tra,e3t,RR_P2,NORM_1,dCa,jpe,
!$omp&        jpf,g0,gdepw,NPPF2,rdt,tra_DIA)      

#ifdef __OPENMP
        mytid = omp_get_thread_num()  ! take the thread ID
#endif
      IF( mytid + jv <= dimen_jvsnu ) THEN

          ji = jarr_snu(1,jv+mytid)
          jj = jarr_snu(2,jv+mytid)

          TOTcalc(jv+mytid) = 0. ! compute total flux due to coccolit.  uptake
          DO jk = 1,jpe          ! loop till the euphotic layer
              TOTcalc(jv+mytid)=TOTcalc(jv+mytid) + RR_P2*NPPF2(ji,jj,jk)*rdt/86400.*e3t(jk)*tmask(ji,jj,jk)
          END DO


C Compute correction factor 
          DO jk=1, jpe
c           dDic= RR_P2*NPPF2(ji,jj,jk)*rdt/86400.*e3t(jk)*tmask(ji,jj,jk)
c           dAlk= 2./12.*RR_P2*NPPF2(ji,jj,jk)*rdt/86400.*e3t(jk)*tmask(ji,jj,jk)
           dDic= RR_P2*NPPF2(ji,jj,jk)*rdt/86400.*tmask(ji,jj,jk)
           dAlk= 2./12.*RR_P2*NPPF2(ji,jj,jk)*rdt/86400.*tmask(ji,jj,jk)
           tra(ji,jj,jk,ppO3c) = tra(ji,jj,jk,ppO3c) - dDic
           tra(ji,jj,jk,ppO3h) = tra(ji,jj,jk,ppO3h) - dAlk
           tra_DIA(ji,jj,jk,ppHT1) = - dDic
           tra_DIA(ji,jj,jk,ppHT2) = - dAlk

          END DO

          g0 =gdepw(jpf)

          DO jk=jpf,jpk
              g1         = gdepw(jk)
              SinkCa(jk) = TOTcalc(jv+mytid)*exp(-(g1-g0)/dCa)
          END DO

          DO jk=jpf+1,jpk
              
              dDic= (SinkCa(jk-1)-SinkCa(jk))/e3t(jk)* tmask(ji,jj,jk)
              dAlk= 2./12.*dDic* tmask(ji,jj,jk)

              tra(ji,jj,jk,ppO3c) = tra(ji,jj,jk,ppO3c) + dDic 
              tra(ji,jj,jk,ppO3h) = tra(ji,jj,jk,ppO3h) + dAlk 
              tra_DIA(ji,jj,jk,ppHT1) = tra_DIA(ji,jj,jk,ppHT1) + dDic
              tra_DIA(ji,jj,jk,ppHT2) = tra_DIA(ji,jj,jk,ppHT2) + dAlk

          END DO


      ENDIF
!$omp end parallel
      END DO ! jv


      END SUBROUTINE hard_tissue_pump
