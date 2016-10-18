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
      integer :: jn_high, jn_on_all
      REAL(8) ::  Miss_val =1.e20
      REAL(8) :: Realcounter, Realcounterp1

! omp variables
      INTEGER :: mytid, ntids

#ifdef __OPENMP1
      INTEGER ::  omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
      EXTERNAL :: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
#endif

      ave_partTime = MPI_WTIME()



#ifdef __OPENMP1
      ntids = omp_get_max_threads() ! take the number of threads
      mytid = -1000000
#else
      ntids = 1
      mytid = 0
#endif


!     FIRST, LOW FREQUENCY
      Realcounter   =    REAL(ave_counter_2  , 8)
      Realcounterp1 = 1./REAL(ave_counter_2+1, 8)

      DO jn=1 ,jptra,ntids

!!!$omp parallel default(none) private(jk,jj,ji,mytid)
!!!$omp&                       shared(jpk,jpj,jpi,jn,tmask,traIO,trn,Miss_val,Realcounter,Realcounterp1)


#ifdef __OPENMP1
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
!!!$omp    end parallel

      END DO




      Realcounter   =    REAL(ave_counter_1  , 8)! ****************** HIGH FREQUENCY
      Realcounterp1 = 1./REAL(ave_counter_1+1, 8)
      DO jn_high=1 ,jptra_high,ntids

!!!$omp parallel default(none) private(jk,jj,ji,mytid,jn_on_all)
!!!$omp& shared(jpk,jpj,jpi,jn_high,jptra_high,highfreq_table,tmask,traIO_HIGH,trn,Miss_val,Realcounter,Realcounterp1)

#ifdef __OPENMP1
           mytid = omp_get_thread_num()  ! take the thread ID
#endif
       IF( mytid + jn_high <= jptra_high ) then
          jn_on_all = highfreq_table(jn_high+mytid)
          DO jk=1, jpk
             DO jj=1, jpj
                DO ji=1, jpi
               IF(tmask(ji,jj,jk) .NE. 0.) THEN
                  traIO_HIGH(ji,jj,jk,jn_high+mytid)= &
     &           (traIO_HIGH(ji,jj,jk,jn_high+mytid)*Realcounter+trn(ji,jj,jk,jn_on_all))*Realcounterp1
               ELSE
                  traIO_HIGH(ji,jj,jk,jn_high+mytid)=Miss_val
               ENDIF
                END DO
             END DO
          END DO
       ENDIF
!!!$omp    end parallel

      END DO


!     *****************  PHYS *****************************************************
      if (freq_ave_phys.eq.1) then
          Realcounter   =    REAL(ave_counter_1  , 8)
          Realcounterp1 = 1./REAL(ave_counter_1+1, 8)
      else
          Realcounter   =    REAL(ave_counter_2  , 8)
          Realcounterp1 = 1./REAL(ave_counter_2+1, 8)
      endif


      DO jk=1, jpk
       DO jj=1, jpj
          DO ji=1, jpi
             IF(tmask(ji,jj,jk) .NE. 0.) THEN
                snIO (ji,jj,jk)=(snIO (ji,jj,jk)*Realcounter+sn (ji,jj,jk))*Realcounterp1
                tnIO (ji,jj,jk)=(tnIO (ji,jj,jk)*Realcounter+tn (ji,jj,jk))*Realcounterp1
                wnIO (ji,jj,jk)=(wnIO (ji,jj,jk)*Realcounter+wn (ji,jj,jk))*Realcounterp1
                avtIO(ji,jj,jk)=(avtIO(ji,jj,jk)*Realcounter+avt(ji,jj,jk))*Realcounterp1
                e3tIO(ji,jj,jk)=(e3tIO(ji,jj,jk)*Realcounter+e3t(ji,jj,jk))*Realcounterp1
             ELSE
                snIO (ji,jj,jk)=Miss_val
                tnIO (ji,jj,jk)=Miss_val
                wnIO (ji,jj,jk)=Miss_val
                avtIO(ji,jj,jk)=Miss_val
                e3tIO(ji,jj,jk)=Miss_val
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
        END DO
      END DO

!     *****************  END PHYS *************************************************


!     *****************  DIAGNOSTICS **********************************************

!     FIRST, LOW FREQUENCY

      Realcounter   =    REAL(ave_counter_2  , 8)
      Realcounterp1 = 1./REAL(ave_counter_2+1, 8)

      DO jn=1, jptra_dia,ntids

!!!$omp parallel default(none) private(jk,jj,ji,mytid)
!!!$omp&   shared(jpk,jpj,jpi,jn,tmask,tra_DIA_IO,tra_DIA,Miss_val,Realcounter,Realcounterp1)

#ifdef __OPENMP1
           mytid = omp_get_thread_num()  ! take the thread ID
#endif
      IF( mytid + jn .LE. jptra_dia ) then

         DO jk=1, jpk
            DO jj=1, jpj
               DO ji=1, jpi
                  IF(tmask(ji,jj,jk) .NE. 0.) THEN
                    tra_DIA_IO(ji,jj,jk,jn+mytid)=(tra_DIA_IO(ji,jj,jk,jn+mytid)*Realcounter+ &
     &              tra_DIA(ji,jj,jk,jn+mytid))*Realcounterp1
                  ELSE
                    tra_DIA_IO(ji,jj,jk,jn+mytid)=Miss_val
                  ENDIF
               END DO
            END DO
         END DO
      ENDIF

!!!$omp    end parallel
      END DO

!     *********************  DIAGNOSTICS 2D **********
      DO jn=1, jptra_dia_2d

            DO jj=1, jpj
               DO ji=1, jpi
                  IF(tmask(ji,jj,1) .NE. 0.) THEN ! Warning ! Tested only for surface
                    tra_DIA_2d_IO(ji,jj,jn)=(tra_DIA_2d_IO(ji,jj,jn)*Realcounter+ &
     &              tra_DIA_2d(ji,jj,jn))*Realcounterp1
                  ELSE
                    tra_DIA_2d_IO(ji,jj,jn)=Miss_val
                  ENDIF
               END DO
            END DO

      END DO




      Realcounter   =    REAL(ave_counter_1  , 8) ! ****************** HIGH FREQUENCY
      Realcounterp1 = 1./REAL(ave_counter_1+1, 8)


      DO jn_high=1, jptra_dia_high,ntids

!!!$omp parallel default(none) private(jk,jj,ji,mytid,jn_on_all)
!!!$omp&   shared(jpk,jpj,jpi,jn_high,jptra_dia_high,highfreq_table_dia, tmask,tra_DIA_IO_HIGH,tra_DIA,Miss_val,
!!!$omp&   Realcounter,Realcounterp1)

#ifdef __OPENMP1
           mytid = omp_get_thread_num()  ! take the thread ID
#endif

      IF (mytid + jn_high .LE. jptra_dia_high)  then
          IF (mytid + jn_high .LE. jptra_dia ) then
             jn_on_all = highfreq_table_dia(jn_high+mytid)

             DO jk=1, jpk
             DO jj=1, jpj
             DO ji=1, jpi
                IF(tmask(ji,jj,jk) .NE. 0.) THEN
                   tra_DIA_IO_HIGH(ji,jj,jk,jn_high+mytid)= &
     &            (tra_DIA_IO_HIGH(ji,jj,jk,jn_high+mytid)*Realcounter+tra_DIA(ji,jj,jk,jn_on_all))*Realcounterp1 
                ELSE
                   tra_DIA_IO_HIGH(ji,jj,jk,jn_high+mytid)=Miss_val
                ENDIF
             END DO
             END DO
             END DO
          ENDIF
      ENDIF
!!!$omp    end parallel
      END DO

!     *********************  DIAGNOSTICS 2D **********

      DO jn_high=1, jptra_dia2d_high
             jn_on_all = highfreq_table_dia2d(jn_high)

             DO jj=1, jpj
             DO ji=1, jpi
                IF(tmask(ji,jj,1) .NE. 0.) THEN
                   tra_DIA_2d_IO_HIGH(ji,jj,jn_high)= &
     &            (tra_DIA_2d_IO_HIGH(ji,jj,jn_high)*Realcounter+tra_DIA_2d(ji,jj,jn_on_all))*Realcounterp1
                ELSE
                   tra_DIA_2d_IO_HIGH(ji,jj,jn_high)=Miss_val
                ENDIF
             END DO
             END DO

      END DO


      ave_partTime = MPI_WTIME() - ave_partTime
      ave_TotTime = ave_TotTime  + ave_partTime

      END SUBROUTINE trcave

