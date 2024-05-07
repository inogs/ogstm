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
      integer :: queue

      ave_partTime = MPI_WTIME()

      queue=1

!     FIRST, LOW FREQUENCY
      elapsed_time         =     elapsed_time_2
      inv_incremented_time = 1./(elapsed_time_2 + rdt)

      !$acc update device(traIO,trn,umask,vmask,tmask,traIO_HIGH,highfreq_table,snIO,tnIO,wnIO,avtIO,e3tIO,unIO,vnIO,sn,tn,wn,avt,e3t,un,vn,tra_DIA_IO,tra_DIA,tra_DIA_2d_IO,tra_DIA_2d,vatmIO,empIO,qsrIO,vatm,emp,qsr,highfreq_table_dia,tra_DIA_IO_HIGH,tra_DIA_2d_IO_HIGH,highfreq_table_dia2d)

      !$acc parallel loop gang vector collapse(4) default(present) async(queue)
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
      !$acc end parallel loop

       



      elapsed_time         =     elapsed_time_1
      inv_incremented_time = 1./(elapsed_time_1 + rdt)! ****************** HIGH FREQUENCY
      !$acc parallel loop gang vector collapse(4) default(present) async(queue)
      DO jn_high=1 ,jptra_high



       
#ifndef _OPENACC
          jn_on_all = highfreq_table(jn_high)
#endif
          DO ji=1, jpi
          DO jj=1, jpj
          DO jk=1, jpk
#ifdef _OPENACC
               jn_on_all = highfreq_table(jn_high)
#endif
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
      !$acc end parallel loop


!     *****************  PHYS *****************************************************
      if (freq_ave_phys.eq.1) then
          elapsed_time         =     elapsed_time_1
          inv_incremented_time = 1./(elapsed_time_1 + rdt)
      else
          elapsed_time         =     elapsed_time_2
          inv_incremented_time = 1./(elapsed_time_2 + rdt)
      endif


          !$acc parallel loop gang vector collapse(3) default(present) async(queue)
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
      !$acc end parallel loop

      !$acc parallel loop gang vector collapse(2) default(present) async(queue)
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
      !$acc end parallel loop


!     *****************  END PHYS *************************************************


!     *****************  DIAGNOSTICS **********************************************

      if (lbfm) THEN
!     FIRST, LOW FREQUENCY

      elapsed_time         =     elapsed_time_2
      inv_incremented_time = 1./(elapsed_time_2 + rdt)


         !$acc parallel loop gang vector collapse(4) default(present) async(queue)
         DO jn = 1,jptra_dia
         DO ji=1, jpi
         DO jj=1, jpj
         DO jk=1, jpk
           IF(tmask(jk,jj,ji) .NE. 0.) THEN
           tra_DIA_IO(jk,jj,ji,jn)=(tra_DIA_IO(jk,jj,ji,jn)*elapsed_time+tra_DIA(jk,jj,ji,jn)*rdt)*inv_incremented_time
           ELSE
                tra_DIA_IO(jk,jj,ji,jn )=Miss_val
           ENDIF
         END DO
         END DO
         END DO
         ENDDO
         !$acc end parallel loop

!     *********************  DIAGNOSTICS 2D **********


          !$acc parallel loop gang vector collapse(2) default(present) async(queue)
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
          !$acc end parallel loop





      elapsed_time         =     elapsed_time_1      ! ****************** HIGH FREQUENCY
      inv_incremented_time = 1./(elapsed_time_1 + rdt)


      if (jptra_dia_high.gt.0) THEN
         !$acc parallel loop gang vector collapse(4) default(present) async(queue)
         DO jn_high=1, jptra_dia_high
#ifndef _OPENACC
         jn_on_all = highfreq_table_dia(jn_high )
#endif
         DO ji=1, jpi
         DO jj=1, jpj
         DO jk=1, jpk
#ifdef _OPENACC
         jn_on_all = highfreq_table_dia(jn_high )
#endif
         IF(tmask(jk,jj,ji) .NE. 0.) THEN
               tra_DIA_IO_HIGH(jk,jj,ji,jn_high )= &
     &         (tra_DIA_IO_HIGH(jk,jj,ji,jn_high )*elapsed_time+tra_DIA(jk,jj,ji,jn_on_all)*rdt)*inv_incremented_time
         ELSE
            tra_DIA_IO_HIGH(jk,jj,ji,jn_high )=Miss_val
         ENDIF
         END DO
         END DO
         END DO
         END DO
         !$acc end parallel loop
      endif


!     *********************  DIAGNOSTICS 2D **********


       if (jptra_dia2d_high.gt.0) THEN
             !$acc parallel loop gang vector collapse(3) default(present) async(queue)
             DO ji=1, jpi
             DO jj=1, jpj
                DO jn_high=1, jptra_dia2d_high
                  IF(tmask(1,jj,ji) .NE. 0.) THEN
                     jn_on_all = highfreq_table_dia2d(jn_high)
                       tra_DIA_2d_IO_HIGH(jn_high,jj,ji)= &
         &            (tra_DIA_2d_IO_HIGH(jn_high,jj,ji)*elapsed_time+tra_DIA_2d(jn_on_all,jj,ji)*rdt)*inv_incremented_time
                  ELSE
                     tra_DIA_2d_IO_HIGH(jn_high,jj,ji)=Miss_val
                  ENDIF
                END DO
             END DO
             END DO
             !$acc end parallel loop
        endif


      endif ! lfbm

      !$acc wait(queue)
      !$acc update host(traIO,traIO_HIGH,snIO,tnIO,wnIO,avtIO,e3tIO,unIO,vnIO,vatmIO,empIO,qsrIO,tra_DIA_IO,tra_DIA_2d_IO,tra_DIA_2d,tra_DIA_IO_HIGH,tra_DIA_2d_IO_HIGH)

      ave_partTime = MPI_WTIME() - ave_partTime
      ave_TotTime = ave_TotTime  + ave_partTime

      END SUBROUTINE trcave

