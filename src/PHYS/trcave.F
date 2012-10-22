      SUBROUTINE trcave
      USE myalloc
      USE IO_mem
      USE FN_mem
#ifdef key_mpp
      USE myalloc_mpp
#endif

      implicit none
!     local
      integer ji,jj,jk,jn
      REAL(8) ::  Miss_val =1.e20
      real(8) Realcounter, Realcounterp1
! omp variables
      INTEGER :: mytid, ntids!, itid

#ifdef __OPENMP
      INTEGER ::  omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
      EXTERNAL :: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
#endif

      ave_partTime = MPI_WTIME()



#ifdef __OPENMP
      ntids = omp_get_max_threads() ! take the number of threads
      mytid = -1000000
#else
      ntids = 1
      mytid = 0
#endif

      Realcounter   =    REAL(ave_counter  , 8)
      Realcounterp1 = 1./REAL(ave_counter+1, 8)




      DO jn=1 ,jptra,ntids
!$omp parallel default(none) private(jk,jj,ji,mytid) shared(jpk,jpj,jpi,jn,tmask,traIO,Realcounter,trn,Realcounterp1,
!$omp&                                                      Miss_val)

#ifdef __OPENMP
           mytid = omp_get_thread_num()  ! take the thread ID
#endif
       IF( mytid + jn <= jptra ) then
          DO jk=1, jpk
             DO jj=1, jpj
                DO ji=1, jpi
               IF(tmask(ji,jj,jk) .NE. 0.) THEN
                  traIO(ji,jj,jk,jn+mytid)=(traIO(ji,jj,jk,jn+mytid)*Realcounter+trn(ji,jj,jk,jn+mytid))*Realcounterp1
               ELSE
                  traIO(ji,jj,jk,jn+mytid)=Miss_val
               ENDIF
                END DO
             END DO
          END DO
       ENDIF
!$omp    end parallel
      END DO

      DO jk=1, jpk
       DO jj=1, jpj
          DO ji=1, jpi
             IF(tmask(ji,jj,jk) .NE. 0.) THEN
                snIO (ji,jj,jk)=(snIO (ji,jj,jk)*Realcounter+sn (ji,jj,jk))*Realcounterp1
                tnIO (ji,jj,jk)=(tnIO (ji,jj,jk)*Realcounter+tn (ji,jj,jk))*Realcounterp1
                wnIO (ji,jj,jk)=(wnIO (ji,jj,jk)*Realcounter+wn (ji,jj,jk))*Realcounterp1
                avtIO(ji,jj,jk)=(avtIO(ji,jj,jk)*Realcounter+avt(ji,jj,jk))*Realcounterp1
             ELSE
                snIO (ji,jj,jk)=Miss_val
                tnIO (ji,jj,jk)=Miss_val
                wnIO (ji,jj,jk)=Miss_val
                avtIO(ji,jj,jk)=Miss_val
             ENDIF


             IF(umask(ji,jj,jk) .NE. 0.) THEN
                unIO(ji,jj,jk)=(unIO(ji,jj,jk)*Realcounter+un(ji,jj,jk))*Realcounterp1
             ELSE
                unIO(ji,jj,jk)=Miss_val
             ENDIF


             IF(vmask(ji,jj,jk) .NE. 0.) THEN
                vnIO(ji,jj,jk)=(vnIO(ji,jj,jk)*Realcounter+vn(ji,jj,jk))*Realcounterp1
             ELSE
                vnIO(ji,jj,jk)=Miss_val
             ENDIF

          END DO
       END DO
      END DO

      DO jj=1, jpj
        DO ji=1, jpi
           IF (tmask(ji,jj,1) .NE. 0.) THEN
               vatmIO(ji,jj)=(vatmIO(ji,jj)*Realcounter+vatm(ji,jj))*Realcounterp1
               empIO (ji,jj)=(empIO (ji,jj)*Realcounter+emp (ji,jj))*Realcounterp1
               qsrIO (ji,jj)=(qsrIO (ji,jj)*Realcounter+qsr (ji,jj))*Realcounterp1
           ELSE
               vatmIO(ji,jj)=Miss_val
               empIO (ji,jj)=Miss_val
               qsrIO (ji,jj)=Miss_val
           ENDIF
           IF(umask(ji,jj,1) .NE. 0.) THEN
               bblxIO(ji,jj)=(bblxIO(ji,jj)*Realcounter+bblx(ji,jj))*Realcounterp1
           ELSE
           bblxIO(ji,jj)=Miss_val
           ENDIF
           IF(vmask(ji,jj,1) .NE. 0.) THEN
               bblyIO(ji,jj)=(bblyIO(ji,jj)*Realcounter+bbly(ji,jj))*Realcounterp1
           ELSE
           bblyIO(ji,jj)=Miss_val
           ENDIF
        END DO
      END DO

      DO jn=1, jptra_dia,ntids
!$omp parallel default(none) private(jk,jj,ji,mytid) shared(jpk,jpj,jpi,jn,tmask,tra_ppIO,Realcounter,tra_pp,
!$omp&   Realcounterp1,Miss_val)

#ifdef __OPENMP
           mytid = omp_get_thread_num()  ! take the thread ID
#endif
      IF( mytid + jn <= jptra_dia ) then
         DO jk=1, jpk
            DO jj=1, jpj
               DO ji=1, jpi
              IF(tmask(ji,jj,jk) .NE. 0.) THEN
                 tra_ppIO(ji,jj,jk,jn+mytid)=(tra_ppIO(ji,jj,jk,jn+mytid)*Realcounter+
     &            tra_pp(ji,jj,jk,jn+mytid))*Realcounterp1
              ELSE
                 tra_ppIO(ji,jj,jk,jn+mytid)=Miss_val
              ENDIF
               END DO
            END DO
         END DO
      ENDIF
!$omp    end parallel
      END DO

! here when nwritetrc  is ok we have to pass the average field!


         ave_counter = ave_counter +1  ! start incrementing our counter

      ave_partTime = MPI_WTIME() - ave_partTime
      ave_TotTime = ave_TotTime  + ave_partTime

      END SUBROUTINE trcave

