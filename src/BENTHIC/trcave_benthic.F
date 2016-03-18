      SUBROUTINE trcave_benthic()
      USE myalloc
      USE BIOBEN_mem
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


!     FIRST, LOW FREQUENCY
      Realcounter   =    REAL(ave_counter_2  , 8)
      Realcounterp1 = 1./REAL(ave_counter_2+1, 8)

      DO jn=1 ,jptra_b,ntids

!$omp parallel default(none) private(jk,jj,ji,mytid)
!$omp&                       shared(jpk_b,jpj,jpi,jn,tmask,traIO,trn,Miss_val,Realcounter,Realcounterp1)


#ifdef __OPENMP
           mytid = omp_get_thread_num()  ! take the thread ID
#endif
       IF( mytid + jn <= jptra_b ) then

          DO jk=1, jpk_ben
             DO jj=1, jpj
                DO ji=1, jpi
               IF(tmask(ji,jj,1) .NE. 0.) THEN
                  benIO(ji,jj,jk,jn+mytid)=(benIO(ji,jj,jk,jn+mytid)*Realcounter+trn_ben(ji,jj,jk,jn+mytid))*Realcounterp1
               ELSE
                  benIO(ji,jj,jk,jn+mytid)=Miss_val
               ENDIF
                END DO
             END DO
          END DO
       ENDIF
!$omp    end parallel

      END DO




           Realcounter   =    REAL(ave_counter_1  , 8)! ****************** HIGH FREQUENCY
           Realcounterp1 = 1./REAL(ave_counter_1+1, 8)
           DO jn_high=1 ,jptra_high_ben,ntids

!$omp parallel default(none) private(jk,jj,ji,mytid,jn_on_all)
!$omp& shared(jpk,jpj,jpi,jn_high,jptra_high,highfreq_table,tmask,traIO_HIGH,trn,Miss_val,Realcounter,Realcounterp1)

#ifdef __OPENMP
           mytid = omp_get_thread_num()  ! take the thread ID
#endif
            IF( mytid + jn_high <= jptra_high_ben ) then
               jn_on_all = highfreq_table_ben(jn_high+mytid)
               DO jk=1, jpk_ben
                  DO jj=1, jpj
                     DO ji=1, jpi
                    IF(tmask(ji,jj,1) .NE. 0.) THEN
                       benIO_HIGH(ji,jj,jk,jn_high+mytid)= 
     &              (benIO_HIGH(ji,jj,jk,jn_high+mytid)*Realcounter+trn_ben(ji,jj,jk,jn_on_all))*Realcounterp1
                    ELSE
                       benIO_HIGH(ji,jj,jk,jn_high+mytid)=Miss_val
                    ENDIF
                     END DO
                  END DO
               END DO
            ENDIF
!$omp    end parallel

           END DO


!     *****************  DIAGNOSTICS **********************************************

!     FIRST, LOW FREQUENCY

           Realcounter   =    REAL(ave_counter_2  , 8)
           Realcounterp1 = 1./REAL(ave_counter_2+1, 8)

           DO jn=1, jptra_dia_b,ntids

!$omp parallel default(none) private(jk,jj,ji,mytid)
!$omp&   shared(jpk_ben,jpj,jpi,jn,tmask,ben_DIA_IO,ben_DIA,Miss_val,Realcounter,Realcounterp1)

#ifdef __OPENMP
           mytid = omp_get_thread_num()  ! take the thread ID
#endif
           IF( mytid + jn .LE. jptra_dia_b ) then

              DO jk=1, jpk_ben
                 DO jj=1, jpj
                    DO ji=1, jpi
                       IF(tmask(ji,jj,1) .NE. 0.) THEN
                         ben_DIA_IO(ji,jj,jk,jn+mytid)=(ben_DIA_IO(ji,jj,jk,jn+mytid)*Realcounter+
     &              ben_DIA(ji,jj,jk,jn+mytid))*Realcounterp1
                       ELSE
                         ben_DIA_IO(ji,jj,jk,jn+mytid)=Miss_val
                       ENDIF
                    END DO
                 END DO
              END DO

           ENDIF

!$omp    end parallel
           END DO

!     *********************  DIAGNOSTICS 2D **********
           DO jn=1, jptra_dia_b_2d

                 DO jj=1, jpj
                    DO ji=1, jpi
                       IF(tmask(ji,jj,1) .NE. 0.) THEN ! Warning ! Tested only for surface
                         ben_DIA_2d_IO(ji,jj,jn)=(ben_DIA_2d_IO(ji,jj,jn)*Realcounter+
     &                   ben_DIA_2d(ji,jj,jn))*Realcounterp1
                       ELSE
                         ben_DIA_2d_IO(ji,jj,jn)=Miss_val
                       ENDIF
                    END DO
                 END DO

           END DO




           Realcounter   =    REAL(ave_counter_1  , 8) ! ****************** HIGH FREQUENCY
           Realcounterp1 = 1./REAL(ave_counter_1+1, 8)


           DO jn_high=1, jptra_dia_high_ben,ntids

!$omp parallel default(none) private(jk,jj,ji,mytid,jn_on_all)
!$omp&   shared(jpk_ben,jpj,jpi,jn_high,jptra_dia_high,highfreq_table_dia, tmask,tra_DIA_IO_HIGH,tra_DIA,Miss_val,
!$omp&   Realcounter,Realcounterp1)

#ifdef __OPENMP
           mytid = omp_get_thread_num()  ! take the thread ID
#endif

           IF (mytid + jn_high .LE. jptra_dia_high_ben)  then
               IF (mytid + jn_high .LE. jptra_dia_b ) then
                  jn_on_all = highfreq_table_dia_ben(jn_high+mytid)

                  DO jk=1, jpk_ben
                  DO jj=1, jpj
                  DO ji=1, jpi
                     IF(tmask(ji,jj,1) .NE. 0.) THEN
                        ben_DIA_IO_HIGH(ji,jj,jk,jn_high+mytid)=
     &                 (ben_DIA_IO_HIGH(ji,jj,jk,jn_high+mytid)*Realcounter+ben_DIA(ji,jj,jk,jn_on_all))*Realcounterp1
                     ELSE
                        ben_DIA_IO_HIGH(ji,jj,jk,jn_high+mytid)=Miss_val
                     ENDIF
                  END DO
                  END DO
                  END DO
               ENDIF
           ENDIF
!$omp    end parallel
           END DO

!     *********************  DIAGNOSTICS 2D **********

           DO jn_high=1, jptra_dia2d_high_ben
                  jn_on_all = highfreq_table_dia2d_ben(jn_high)

                  DO jj=1, jpj
                  DO ji=1, jpi
                     IF(tmask(ji,jj,1) .NE. 0.) THEN
                        ben_DIA_2d_IO_HIGH(ji,jj,jn_high)=
     &                 (ben_DIA_2d_IO_HIGH(ji,jj,jn_high)*Realcounter+ben_DIA_2d(ji,jj,jn_on_all))*Realcounterp1
                     ELSE
                        ben_DIA_2d_IO_HIGH(ji,jj,jn_high)=Miss_val
                     ENDIF
                  END DO
                  END DO

           END DO


      ave_partTime = MPI_WTIME() - ave_partTime
      ave_TotTime = ave_TotTime  + ave_partTime

      END SUBROUTINE trcave_benthic

