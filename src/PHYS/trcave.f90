      SUBROUTINE trcave
      USE myalloc
      USE IO_mem
      USE FN_mem
      use mpi

      implicit none
!     local
      integer jk,jj,ji,jn
      integer :: jn_high, jn_on_all
      double precision ::  Miss_val =1.e20
      double precision :: Realcounter, Realcounterp1

      ave_partTime = MPI_WTIME()


!     FIRST, LOW FREQUENCY
      Realcounter   =    REAL(ave_counter_2  , 8)
      Realcounterp1 = 1./REAL(ave_counter_2+1, 8)

      DO jn=1 ,jptra

!!!$omp parallel default(none) private(jk,jj,ji, )
!!!$omp&                       shared(jpk,jpj,jpi,jn,tmask,traIO,trn,Miss_val,Realcounter,Realcounterp1)


                DO ji=1, jpi
             DO jj=1, jpj
          DO jk=1, jpk
               IF(tmask(jk,jj,ji) .NE. 0.) THEN
                  traIO(jk,jj,ji,jn )=(traIO(jk,jj,ji,jn )*Realcounter+trn(jk,jj,ji,jn ))*Realcounterp1
               ELSE
                  traIO(jk,jj,ji,jn )=Miss_val
               ENDIF
                END DO
             END DO
          END DO

      END DO

      ! OPEN(UNIT=10009, FILE='s9.txt', FORM='FORMATTED')
      ! DO jn=1,jptra; DO jk = 1,jpk; DO jj = 1,jpj; DO ji = 1,jpi;
      ! WRITE(10009,200),'S9',jn,jk,jj,ji,traIO(jk,jj,ji,jn)
      ! ENDDO;ENDDO;ENDDO;ENDDO;CLOSE(10009)

       



      Realcounter   =    REAL(ave_counter_1  , 8)! ****************** HIGH FREQUENCY
      Realcounterp1 = 1./REAL(ave_counter_1+1, 8)
      DO jn_high=1 ,jptra_high

!!!$omp parallel default(none) private(jk,jj,ji, ,jn_on_all)
!!!$omp& shared(jpk,jpj,jpi,jn_high,jptra_high,highfreq_table,tmask,traIO_HIGH,trn,Miss_val,Realcounter,Realcounterp1)


       
          jn_on_all = highfreq_table(jn_high)
                DO ji=1, jpi
             DO jj=1, jpj
          DO jk=1, jpk
               IF(tmask(jk,jj,ji) .NE. 0.) THEN
                  traIO_HIGH(jk,jj,ji,jn_high  )= &
     &           (traIO_HIGH(jk,jj,ji,jn_high  )*Realcounter+trn(jk,jj,ji,jn_on_all))*Realcounterp1
               ELSE
                  traIO_HIGH(jk,jj,ji,jn_high  )=Miss_val
               ENDIF
                END DO
             END DO
          END DO
   
!!!$omp    end parallel

      END DO

      ! OPEN(UNIT=10010, FILE='s10.txt', FORM='FORMATTED')
      ! DO jn=1,jptra; DO jk = 1,jpk; DO jj = 1,jpj; DO ji = 1,jpi;
      ! WRITE(10010,200),'S10',jn,jk,jj,ji,traIO(jk,jj,ji,jn)
      ! ENDDO;ENDDO;ENDDO;ENDDO;CLOSE(10010)

!     *****************  PHYS *****************************************************
      if (freq_ave_phys.eq.1) then
          Realcounter   =    REAL(ave_counter_1  , 8)
          Realcounterp1 = 1./REAL(ave_counter_1+1, 8)
      else
          Realcounter   =    REAL(ave_counter_2  , 8)
          Realcounterp1 = 1./REAL(ave_counter_2+1, 8)
      endif


          DO ji=1, jpi
       DO jj=1, jpj
      DO jk=1, jpk
             IF(tmask(jk,jj,ji) .NE. 0.) THEN
                snIO (jk,jj,ji)=(snIO (jk,jj,ji)*Realcounter+sn (jk,jj,ji))*Realcounterp1
                tnIO (jk,jj,ji)=(tnIO (jk,jj,ji)*Realcounter+tn (jk,jj,ji))*Realcounterp1
                wnIO (jk,jj,ji)=(wnIO (jk,jj,ji)*Realcounter+wn (jk,jj,ji))*Realcounterp1
                avtIO(jk,jj,ji)=(avtIO(jk,jj,ji)*Realcounter+avt(jk,jj,ji))*Realcounterp1
                e3tIO(jk,jj,ji)=(e3tIO(jk,jj,ji)*Realcounter+e3t(jk,jj,ji))*Realcounterp1
             ELSE
                snIO (jk,jj,ji)=Miss_val
                tnIO (jk,jj,ji)=Miss_val
                wnIO (jk,jj,ji)=Miss_val
                avtIO(jk,jj,ji)=Miss_val
                e3tIO(jk,jj,ji)=Miss_val
             ENDIF


             IF(umask(jk,jj,ji) .NE. 0.) THEN
                unIO(jk,jj,ji)=(unIO(jk,jj,ji)*Realcounter+un(jk,jj,ji))*Realcounterp1
             ELSE
                unIO(jk,jj,ji)=Miss_val
             ENDIF


             IF(vmask(jk,jj,ji) .NE. 0.) THEN
                vnIO(jk,jj,ji)=(vnIO(jk,jj,ji)*Realcounter+vn(jk,jj,ji))*Realcounterp1
             ELSE
                vnIO(jk,jj,ji)=Miss_val
             ENDIF

          END DO
       END DO
      END DO

      DO jj=1, jpj
        DO ji=1, jpi
           IF (tmask(1,jj,ji) .NE. 0.) THEN
               vatmIO(jj,ji)=(vatmIO(jj,ji)*Realcounter+vatm(jj,ji))*Realcounterp1
       !        write(*,200),ji,jj,vatm(jj,ji),vatmIO(jj,ji)
               empIO (jj,ji)=(empIO (jj,ji)*Realcounter+emp (jj,ji))*Realcounterp1
               qsrIO (jj,ji)=(qsrIO (jj,ji)*Realcounter+qsr (jj,ji))*Realcounterp1
           ELSE
               vatmIO(jj,ji)=Miss_val
               empIO (jj,ji)=Miss_val
               qsrIO (jj,ji)=Miss_val
           ENDIF
        END DO
      END DO
      !STOP


!     *****************  END PHYS *************************************************


!     *****************  DIAGNOSTICS **********************************************

!     FIRST, LOW FREQUENCY

      Realcounter   =    REAL(ave_counter_2  , 8)
      Realcounterp1 = 1./REAL(ave_counter_2+1, 8)

      DO jn=1, jptra_dia

!!!$omp parallel default(none) private(jk,jj,ji, )
!!!$omp&   shared(jpk,jpj,jpi,jn,tmask,tra_DIA_IO,tra_DIA,Miss_val,Realcounter,Realcounterp1)

      !IF( jn .LE. jptra_dia ) then

               DO ji=1, jpi
            DO jj=1, jpj
         DO jk=1, jpk
                  IF(tmask(jk,jj,ji) .NE. 0.) THEN
                    tra_DIA_IO(jn,jk,jj,ji )=(tra_DIA_IO(jn,jk,jj,ji )*Realcounter+ &
     &              tra_DIA(jn,jk,jj,ji))*Realcounterp1
                  ELSE
                    tra_DIA_IO(jn,jk,jj,ji )=Miss_val
                  ENDIF
               END DO
            END DO
         END DO
      !ENDIF
      

!!!$omp    end parallel
      END DO
      ! print *,"---------2",tra_DIA_IO(1,30,15,:)
      ! OPEN(UNIT=10011, FILE='s11.txt', FORM='FORMATTED')
      ! DO jn=1,jptra_dia; DO jk = 1,jpk; DO jj = 1,jpj; DO ji = 1,jpi;
      ! WRITE(10011,200),'S11',jn,jk,jj,ji,tra_DIA_IO(jn,jk,jj,ji)
      ! ENDDO;ENDDO;ENDDO;ENDDO;CLOSE(10011)

!     *********************  DIAGNOSTICS 2D **********
      DO jn=1, jptra_dia_2d

               DO ji=1, jpi
            DO jj=1, jpj
                  IF(tmask(1,jj,ji) .NE. 0.) THEN ! Warning ! Tested only for surface
                    tra_DIA_2d_IO(jj,ji,jn)=(tra_DIA_2d_IO(jj,ji,jn)*Realcounter+ &
     &              tra_DIA_2d(jj,ji,jn))*Realcounterp1
                  ELSE
                    tra_DIA_2d_IO(jj,ji,jn)=Miss_val
                  ENDIF
               END DO
            END DO

      END DO





      Realcounter   =    REAL(ave_counter_1  , 8) ! ****************** HIGH FREQUENCY
      Realcounterp1 = 1./REAL(ave_counter_1+1, 8)


      DO jn_high=1, jptra_dia_high

!!!$omp parallel default(none) private(jk,jj,ji, ,jn_on_all)
!!!$omp&   shared(jpk,jpj,jpi,jn_high,jptra_dia_high,highfreq_table_dia, tmask,tra_DIA_IO_HIGH,tra_DIA,Miss_val,
!!!$omp&   Realcounter,Realcounterp1)



     
          IF (jn_high .LE. jptra_dia_high)  then
          IF (jn_high .LE. jptra_dia ) then
             jn_on_all = highfreq_table_dia(jn_high )

             DO ji=1, jpi
             DO jj=1, jpj
             DO jk=1, jpk
                IF(tmask(jk,jj,ji) .NE. 0.) THEN
                   tra_DIA_IO_HIGH(jk,jj,ji,jn_high )= &
     &            (tra_DIA_IO_HIGH(jk,jj,ji,jn_high )*Realcounter+tra_DIA(jn_on_all,jk,jj,ji))*Realcounterp1 
                ELSE
                   tra_DIA_IO_HIGH(jk,jj,ji,jn_high )=Miss_val
                ENDIF
             END DO
             END DO
             END DO
          ENDIF
          ENDIF
     
!!!$omp    end parallel
      END DO

      ! OPEN(UNIT=10012, FILE='s12.txt', FORM='FORMATTED')
      ! DO jn=1,jptra_dia_high; DO jk = 1,jpk; DO jj = 1,jpj; DO ji = 1,jpi;
      ! WRITE(10012,200),'S12',jn,jk,jj,ji,tra_DIA_IO_HIGH(jk,jj,ji,jn)
      ! ENDDO;ENDDO;ENDDO;ENDDO;CLOSE(10012)

!     *********************  DIAGNOSTICS 2D **********

      DO jn_high=1, jptra_dia2d_high
             jn_on_all = highfreq_table_dia2d(jn_high)

             DO ji=1, jpi
             DO jj=1, jpj
                IF(tmask(1,jj,ji) .NE. 0.) THEN
                   tra_DIA_2d_IO_HIGH(jj,ji,jn_high)= &
     &            (tra_DIA_2d_IO_HIGH(jj,ji,jn_high)*Realcounter+tra_DIA_2d(jj,ji,jn_on_all))*Realcounterp1
                ELSE
                   tra_DIA_2d_IO_HIGH(jj,ji,jn_high)=Miss_val
                ENDIF
             END DO
             END DO

      END DO
      

      ave_partTime = MPI_WTIME() - ave_partTime
      ave_TotTime = ave_TotTime  + ave_partTime

      END SUBROUTINE trcave

