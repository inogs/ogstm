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
      double precision :: elapsed_time, inv_incremented_time

      ave_partTime = MPI_WTIME()


!     FIRST, LOW FREQUENCY
      elapsed_time         =     elapsed_time_2
      inv_incremented_time = 1./(elapsed_time_2 + rdt)

      DO jn=1 ,jptra

          DO ji=1, jpi
          DO jj=1, jpj
          DO jk=1, jpk
               IF(tmask(jk,jj,ji) .NE. 0.) THEN
                  traIO(jk,jj,ji,jn ) = (traIO(jk,jj,ji,jn)*elapsed_time + trn(jk,jj,ji,jn)*rdt)*inv_incremented_time
               ELSE
                  traIO(jk,jj,ji,jn )=Miss_val
               ENDIF
          END DO
          END DO
          END DO

      END DO

       



      elapsed_time         =     elapsed_time_1
      inv_incremented_time = 1./(elapsed_time_1 + rdt)! ****************** HIGH FREQUENCY
      DO jn_high=1 ,jptra_high



       
          jn_on_all = highfreq_table(jn_high)
          DO ji=1, jpi
          DO jj=1, jpj
          DO jk=1, jpk
               IF(tmask(jk,jj,ji) .NE. 0.) THEN
                  traIO_HIGH(jk,jj,ji,jn_high  )= &
     &           (traIO_HIGH(jk,jj,ji,jn_high )*elapsed_time+trn(jk,jj,ji,jn_on_all)*rdt)*inv_incremented_time
               ELSE
                  traIO_HIGH(jk,jj,ji,jn_high  )=Miss_val
               ENDIF
          END DO
          END DO
          END DO
   

      END DO


!     *****************  PHYS *****************************************************
      if (freq_ave_phys.eq.1) then
          elapsed_time         =     elapsed_time_1
          inv_incremented_time = 1./(elapsed_time_1 + rdt)
      else
          elapsed_time         =     elapsed_time_2
          inv_incremented_time = 1./(elapsed_time_2 + rdt)
      endif


          DO ji=1, jpi
       DO jj=1, jpj
      DO jk=1, jpk
             IF(tmask(jk,jj,ji) .NE. 0.) THEN
                snIO (jk,jj,ji)=(snIO (jk,jj,ji)*elapsed_time+sn (jk,jj,ji)*rdt)*inv_incremented_time
                tnIO (jk,jj,ji)=(tnIO (jk,jj,ji)*elapsed_time+tn (jk,jj,ji)*rdt)*inv_incremented_time
                wnIO (jk,jj,ji)=(wnIO (jk,jj,ji)*elapsed_time+wn (jk,jj,ji)*rdt)*inv_incremented_time
                avtIO(jk,jj,ji)=(avtIO(jk,jj,ji)*elapsed_time+avt(jk,jj,ji)*rdt)*inv_incremented_time
                e3tIO(jk,jj,ji)=(e3tIO(jk,jj,ji)*elapsed_time+e3t(jk,jj,ji)*rdt)*inv_incremented_time
             ELSE
                snIO (jk,jj,ji)=Miss_val
                tnIO (jk,jj,ji)=Miss_val
                wnIO (jk,jj,ji)=Miss_val
                avtIO(jk,jj,ji)=Miss_val
                e3tIO(jk,jj,ji)=Miss_val
             ENDIF


             IF(umask(jk,jj,ji) .NE. 0.) THEN
                unIO(jk,jj,ji)=(unIO(jk,jj,ji)*elapsed_time+un(jk,jj,ji)*rdt)*inv_incremented_time
             ELSE
                unIO(jk,jj,ji)=Miss_val
             ENDIF


             IF(vmask(jk,jj,ji) .NE. 0.) THEN
                vnIO(jk,jj,ji)=(vnIO(jk,jj,ji)*elapsed_time+vn(jk,jj,ji)*rdt)*inv_incremented_time
             ELSE
                vnIO(jk,jj,ji)=Miss_val
             ENDIF

          END DO
       END DO
      END DO

      DO jj=1, jpj
        DO ji=1, jpi
           IF (tmask(1,jj,ji) .NE. 0.) THEN
               vatmIO(jj,ji)=(vatmIO(jj,ji)*elapsed_time+vatm(jj,ji)*rdt)*inv_incremented_time
               empIO (jj,ji)=(empIO (jj,ji)*elapsed_time+emp (jj,ji)*rdt)*inv_incremented_time
               qsrIO (jj,ji)=(qsrIO (jj,ji)*elapsed_time+qsr (jj,ji)*rdt)*inv_incremented_time
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

      elapsed_time         =     elapsed_time_2
      inv_incremented_time = 1./(elapsed_time_2 + rdt)



         DO ji=1, jpi
         DO jj=1, jpj
         DO jk=1, jpk
           IF(tmask(jk,jj,ji) .NE. 0.) THEN
           tra_DIA_IO(:,jk,jj,ji)=(tra_DIA_IO(:,jk,jj,ji)*elapsed_time+tra_DIA(:,jk,jj,ji)*rdt)*inv_incremented_time
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
                    tra_DIA_2d_IO(:,jj,ji)=(tra_DIA_2d_IO(:,jj,ji)*elapsed_time+ &
     &              tra_DIA_2d(:,jj,ji)*rdt)*inv_incremented_time
                  ELSE
                    tra_DIA_2d_IO(:,jj,ji)=Miss_val
                  ENDIF
          END DO
          END DO





      elapsed_time         =     elapsed_time_1      ! ****************** HIGH FREQUENCY
      inv_incremented_time = 1./(elapsed_time_1 + rdt)


     if (jptra_dia_high.gt.0) THEN
         DO ji=1, jpi
         DO jj=1, jpj
         DO jk=1, jpk
         IF(tmask(jk,jj,ji) .NE. 0.) THEN
            DO jn_high=1, jptra_dia_high
               jn_on_all = highfreq_table_dia(jn_high )
               tra_DIA_IO_HIGH(jn_high, jk,jj,ji )= &
     &         (tra_DIA_IO_HIGH(jn_high, jk,jj,ji )*elapsed_time+tra_DIA(jn_on_all,jk,jj,ji)*rdt)*inv_incremented_time
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
     &            (tra_DIA_2d_IO_HIGH(jn_high,jj,ji)*elapsed_time+tra_DIA_2d(jn_on_all,jj,ji)*rdt)*inv_incremented_time
                END DO
                ELSE
                   tra_DIA_2d_IO_HIGH(:,jj,ji)=Miss_val
                ENDIF
             END DO
             END DO
        endif


      endif ! lfbm

      ave_partTime = MPI_WTIME() - ave_partTime
      ave_TotTime = ave_TotTime  + ave_partTime

      END SUBROUTINE trcave

