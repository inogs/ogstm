      SUBROUTINE trcave
      USE myalloc
      USE IO_mem
      USE FN_mem
      USE OPT_mem
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
      Realcounterp1 = 1.0d0/REAL(ave_counter_2+1, 8)

      DO jn=1 ,jptra



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

       



      Realcounter   =    REAL(ave_counter_1  , 8)! ****************** HIGH FREQUENCY
      Realcounterp1 = 1.0d0/REAL(ave_counter_1+1, 8)
      DO jn_high=1 ,jptra_high



       
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
   

      END DO


!     *****************  PHYS *****************************************************
      if (freq_ave_phys.eq.1) then
          Realcounter   =    REAL(ave_counter_1  , 8)
          Realcounterp1 = 1.0d0/REAL(ave_counter_1+1, 8)
      else
          Realcounter   =    REAL(ave_counter_2  , 8)
          Realcounterp1 = 1.0d0/REAL(ave_counter_2+1, 8)
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


!     *****************  END PHYS *************************************************


!     *****************  DIAGNOSTICS **********************************************

      if (lbfm) THEN
!     FIRST, LOW FREQUENCY

      Realcounter   =    REAL(ave_counter_2  , 8)
      Realcounterp1 = 1.0d0/REAL(ave_counter_2+1, 8)




               DO ji=1, jpi
            DO jj=1, jpj
         DO jk=1, jpk
                  IF(tmask(jk,jj,ji) .NE. 0.) THEN
                    tra_DIA_IO(:,jk,jj,ji ) = (tra_DIA_IO(:,jk,jj,ji )*Realcounter+ &
     &              tra_DIA(:,jk,jj,ji))*Realcounterp1
                  ELSE
                    tra_DIA_IO(:,jk,jj,ji )=Miss_val
                  ENDIF
               END DO
            END DO
      

      END DO

!     *********************  DIAGNOSTICS 2D **********


               DO ji=1, jpi
            DO jj=1, jpj
                  IF(tmask(1,jj,ji) .NE. 0.) THEN ! Warning ! Tested only for surface
                    tra_DIA_2d_IO(:,jj,ji)=(tra_DIA_2d_IO(:,jj,ji)*Realcounter+ &
     &              tra_DIA_2d(:,jj,ji))*Realcounterp1
                  ELSE
                    tra_DIA_2d_IO(:,jj,ji)=Miss_val
                  ENDIF
               END DO
            END DO






      Realcounter   =    REAL(ave_counter_1  , 8) ! ****************** HIGH FREQUENCY
      Realcounterp1 = 1.0d0/REAL(ave_counter_1+1, 8)




     if (jptra_dia_high.gt.0) THEN
             DO ji=1, jpi
             DO jj=1, jpj
             DO jk=1, jpk
                IF(tmask(jk,jj,ji) .NE. 0.) THEN
                DO jn_high=1, jptra_dia_high
                   jn_on_all = highfreq_table_dia(jn_high )
                   tra_DIA_IO_HIGH(jn_high, jk,jj,ji )= &
     &            (tra_DIA_IO_HIGH(jn_high, jk,jj,ji )*Realcounter+tra_DIA(jn_on_all,jk,jj,ji))*Realcounterp1
                 END DO
                ELSE
                   tra_DIA_IO_HIGH(:, jk,jj,ji )=Miss_val
                ENDIF
             END DO
             END DO
             END DO
     endif


!     *********************  DIAGNOSTICS 2D **********


       if (jptra_dia2d_high.gt.0) THEN
             DO ji=1, jpi
             DO jj=1, jpj
                IF(tmask(1,jj,ji) .NE. 0.) THEN
                DO jn_high=1, jptra_dia2d_high
                   jn_on_all = highfreq_table_dia2d(jn_high)
                   tra_DIA_2d_IO_HIGH(jn_high,jj,ji)= &
     &            (tra_DIA_2d_IO_HIGH(jn_high,jj,ji)*Realcounter+tra_DIA_2d(jn_on_all,jj,ji))*Realcounterp1
                END DO
                ELSE
                   tra_DIA_2d_IO_HIGH(:,jj,ji)=Miss_val
                ENDIF
             END DO
             END DO
        endif


      endif ! lfbm
!     *********************  DIAGNOSTICS RT **********

             DO ji=1, jpi
             DO jj=1, jpj
             DO jk=1, jpk
                IF(tmask(jk,jj,ji) .NE. 0.) THEN
                DO jn_high=1, nlt
                   Ed_DIA_IO(jk,jj,ji,jn_high )= &
     &            (Ed_DIA_IO(jk,jj,ji,jn_high)*Realcounter+Ed(jk,jj,ji,jn_high))*Realcounterp1

                   Es_DIA_IO(jk,jj,ji,jn_high )= &
     &            (Es_DIA_IO(jk,jj,ji,jn_high)*Realcounter+Es(jk,jj,ji,jn_high))*Realcounterp1

                   Eu_DIA_IO(jk,jj,ji,jn_high )= &
     &            (Eu_DIA_IO(jk,jj,ji,jn_high)*Realcounter+Eu(jk,jj,ji,jn_high))*Realcounterp1
                 END DO
                ELSE
                   Ed_DIA_IO(jk,jj,ji,: )=Miss_val
                   Es_DIA_IO(jk,jj,ji,: )=Miss_val
                   Eu_DIA_IO(jk,jj,ji,: )=Miss_val
                ENDIF
             END DO
             END DO
             END DO

      ave_partTime = MPI_WTIME() - ave_partTime
      ave_TotTime = ave_TotTime  + ave_partTime

      END SUBROUTINE trcave

