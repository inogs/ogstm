
      SUBROUTINE hard_tissue_pump()
#ifdef BFMv2
!!!---------------------------------------------------------------------
!!!
!!!                       ROUTINE hard_tissue_pump
!!!                     *******************
!!!
!!!  PURPOSE :
!!!  ---------
!!!	compensate along the water column to preserve total mass
!!!
!!!
!!   METHOD :
!!   -------
!!    Must be called between trczdf and trcnxt

!!----------------------------------------------------------------------

       USE myalloc
       USE FN_mem
       USE BIO_mem
       
       IMPLICIT NONE


!!----------------------------------------------------------------------
!! local declarations
!! ==================

! omp variables
      INTEGER :: mytid, ntids, itid

#ifdef __OPENMP
      INTEGER ::  omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
      EXTERNAL :: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
#endif

      INTEGER jk,jj,ji,jn,jv
      INTEGER jpe,jpf
      REAL(8) RR_P2
      REAL(8) NORM_1, dCa
      REAL(8) g0,g1,dDic,dAlk
      REAL(8) SinkCa(jpk)


!!----------------------------------------------------------------------
!! statement functions
!! ===================


#ifdef __OPENMP
      ntids = omp_get_max_threads() ! take the number of threads
      mytid = -1000000
#else
      ntids = 1
      mytid = 0
#endif

! ***************************************************
      jpe    = jpk_eu
      jpf    = 2 ! CaCO3 dissolution along the whole water column  
      RR_P2  = 0.2   ![]
      dCa    = 1500. ![m]   Mediterranean Sea mean depth
      SinkCa = 0.
! Compute vertical integral
      DO jv  = 1,dimen_jvsnu, ntids ! cicla sui punti terra 2D

!$omp   parallel default(none) private(mytid,jk,jj,ji,jn,g1,SinkCa,dDic,dAlk)
!$omp&  shared(dimen_jvsnu,jarr_snu,jv,jpk,tmask,TOTcalc,tra,e3t,RR_P2,NORM_1,dCa,jpe,
!$omp&        jpf,g0,gdepw,NPPF2,rdt,tra_DIA)      

#ifdef __OPENMP
        mytid = omp_get_thread_num()  ! take the thread ID
#endif
      IF( mytid + jv <= dimen_jvsnu ) THEN

          jj = jarr_snu(1,jv+mytid)
          ji = jarr_snu(2,jv+mytid)

          TOTcalc(jv+mytid) = 0. ! compute total flux due to coccolit.  uptake
          DO jk = 1,jpe          ! loop till the euphotic layer
              TOTcalc(jv+mytid)=TOTcalc(jv+mytid) + RR_P2*NPPF2(jk,jj,ji)*rdt/86400.*e3t(jk,jj,ji)*tmask(jk,jj,ji)
          END DO


! Compute correction factor 
          DO jk=1, jpe
           dDic= RR_P2*NPPF2(jk,jj,ji)*rdt/86400.*tmask(jk,jj,ji)
           dAlk= 2./12.*RR_P2*NPPF2(jk,jj,ji)*rdt/86400.*tmask(jk,jj,ji)
           tra(jk,jj,ji,ppO3c) = tra(jk,jj,ji,ppO3c) - dDic
           tra(jk,jj,ji,ppO3h) = tra(jk,jj,ji,ppO3h) - dAlk
           tra_DIA(ppHT1,jk,jj,ji) = - dDic
           tra_DIA(ppHT2,jk,jj,ji) = - dAlk


          END DO

          g0 =gdepw(jpf)

          DO jk=jpf,jpk
              g1         = gdepw(jk)
              SinkCa(jk) = TOTcalc(jv+mytid)*exp(-(g1-g0)/dCa)
          END DO

          DO jk=jpf+1,jpk
              
              dDic= (SinkCa(jk-1)-SinkCa(jk))/e3t(jk,jj,ji)* tmask(jk,jj,ji)
              dAlk= 2./12.*dDic* tmask(jk,jj,ji)

              tra(jk,jj,ji,ppO3c) = tra(jk,jj,ji,ppO3c) + dDic 
              tra(jk,jj,ji,ppO3h) = tra(jk,jj,ji,ppO3h) + dAlk 
              tra_DIA(ppHT1,jk,jj,ji) = tra_DIA(ppHT1,jk,jj,ji) + dDic
              tra_DIA(ppHT2,jk,jj,ji) = tra_DIA(ppHT2,jk,jj,ji) + dAlk

          END DO


      ENDIF
!$omp end parallel
      END DO ! jv

#endif
      END SUBROUTINE hard_tissue_pump
