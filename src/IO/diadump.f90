      SUBROUTINE diadump(datemean,datefrom,dateTo,FREQ_GROUP)
!     ******************
      USE calendar
      USE myalloc
      USE IO_mem
      USE FN_mem
      USE TIME_MANAGER
      use mpi
      USE ogstm_mpi_module
      USE MPI_GATHER_INFO
      USE MATRIX_VARS
      USE NODES_MODULE
      USE DTYPE_PROCS_STRING_MODULE

      IMPLICIT NONE


      CHARACTER(LEN=17), INTENT(IN) :: datemean, dateFrom, dateTo
      INTEGER, INTENT(IN) :: FREQ_GROUP
      INTEGER jk,jj,ji
      INTEGER ind, i_contribution, j_contribution
      CHARACTER(10) newfile
      CHARACTER(LEN=42) forcing_file
      CHARACTER(LEN=60) bkpname
      CHARACTER(LEN=11) DIR
      logical IsBackup
      double precision :: elapsed_time


      CHARACTER(LEN=56) dia_file_nc
      CHARACTER(LEN=20)  var

      INTEGER idrank, ierr, istart, jstart, iPe, iPd, jPe, jPd, status(MPI_STATUS_SIZE)
      INTEGER irange, jrange
      INTEGER totistart, totiend, relistart, reliend
      INTEGER totjstart, totjend, reljstart, reljend
      double precision ::  Miss_val =1.e20
      INTEGER :: nVars, counter_var_2d, counter_var_high_2d,counter_var_diag, counter_var_diag_high
      CHARACTER(LEN=20) ::  var_to_store_diag_2d, var_to_store_diag
      INTEGER :: n_dumping_cycles, jv, ivar, writing_rank, ind_col
      INTEGER :: var_to_send_2D, var_high_to_send_2D
      INTEGER :: var_to_send, var_high_to_send
     ! call mppsync()
! ----------------------------------------
      IsBackup =  (datemean.eq.dateTo)
      if (lwp) write(*,*) 'diadump IsBackup = ',IsBackup, ' group ' ,FREQ_GROUP
! ----------------------------------------
      bkpname  = 'ave.20111231-15:30:00.N1p.nc.bkp'
      if (IsBackup) then
         forcing_file   = 'AVE_PHYS/ave.'//datemean//'.phys.nc.bkp'
      else
         forcing_file   = 'AVE_PHYS/ave.'//datemean//'.phys.nc'
      endif

        write(*,*) 'old phys'

      SELECT CASE (FREQ_GROUP)
        CASE (1) 
       elapsed_time=elapsed_time_1
       DIR='AVE_FREQ_1/'
        CASE (2) 
       elapsed_time=elapsed_time_2
       DIR='AVE_FREQ_2/'
      END SELECT

      if (lwp) then
          totsnIO   = Miss_val
          tottnIO   = Miss_val
          totunIO   = Miss_val
          totvnIO   = Miss_val
          totwnIO   = Miss_val
          totavtIO  = Miss_val
          tote3tIO  = Miss_val
          totvatmIO = Miss_val
          totempIO  = Miss_val
          totqsrIO  = Miss_val
      endif

!      PHYSICS FIRST!!
      if ( freq_ave_phys.eq.FREQ_GROUP) then
      ! *************** START COLLECTING DATA *****************************
       if (myrank == 0) then                    ! IF LABEL 1


! ******* myrank 0 sets indexes of tot matrix where to place its own part

          iPd    = nldi
          iPe    = nlei
          jPd    = nldj
          jPe    = nlej
          istart = nimpp
          jstart = njmpp
          irange    = iPe - iPd + 1
          jrange    = jPe - jPd + 1
          totistart = istart + iPd - 1 
          totiend   = totistart + irange - 1
          totjstart = jstart + jPd - 1 
          totjend   = totjstart + jrange - 1
          relistart = 1 + iPd - 1      
          reliend   = relistart + irange - 1
          reljstart = 1 + jPd - 1      
          reljend   = reljstart + jrange - 1

          totsnIO  (:, totjstart:totjend,totistart:totiend) = snIO    (:, reljstart:reljend,relistart:reliend)
          tottnIO  (:, totjstart:totjend,totistart:totiend) = tnIO    (:, reljstart:reljend,relistart:reliend)
          totvatmIO(totjstart:totjend,totistart:totiend) = vatmIO  (reljstart:reljend,relistart:reliend)
          totempIO (totjstart:totjend,totistart:totiend) = empIO   (reljstart:reljend,relistart:reliend)
          totqsrIO (totjstart:totjend,totistart:totiend) = qsrIO   (reljstart:reljend,relistart:reliend)
          totunIO  (:, totjstart:totjend,totistart:totiend) = unIO    (:, reljstart:reljend,relistart:reliend)
          totvnIO  (:, totjstart:totjend,totistart:totiend) = vnIO    (:, reljstart:reljend,relistart:reliend)
          totwnIO  (:, totjstart:totjend,totistart:totiend) = wnIO    (:, reljstart:reljend,relistart:reliend)
          totavtIO (:, totjstart:totjend,totistart:totiend) = avtIO   (:, reljstart:reljend,relistart:reliend)
          tote3tIO (:, totjstart:totjend,totistart:totiend) = e3tIO   (:, reljstart:reljend,relistart:reliend)

           do idrank = 1,mpi_glcomm_size-1
 ! **************  myrank 0 is receiving from the others their buffer  ****
               call MPI_RECV(jpi_rec    , 1,                 mpi_integer, idrank, 1,mpi_comm_world, status, ierr) !* first info to know where idrank is working
               call MPI_RECV(jpj_rec    , 1,                 mpi_integer, idrank, 2,mpi_comm_world, status, ierr)
               call MPI_RECV(istart     , 1,                 mpi_integer, idrank, 3,mpi_comm_world, status, ierr)
               call MPI_RECV(jstart     , 1,                 mpi_integer, idrank, 4,mpi_comm_world, status, ierr)
               call MPI_RECV(iPe        , 1,                 mpi_integer, idrank, 5,mpi_comm_world, status, ierr)
               call MPI_RECV(jPe        , 1,                 mpi_integer, idrank, 6,mpi_comm_world, status, ierr)
               call MPI_RECV(iPd        , 1,                 mpi_integer, idrank, 7,mpi_comm_world, status, ierr)
               call MPI_RECV(jPd        , 1                 ,mpi_integer, idrank, 8,mpi_comm_world, status, ierr)
       call MPI_RECV(buffsn  ,jpi_rec*jpj_rec*jpk          ,mpi_real8,idrank, 11,mpi_comm_world, status, ierr)
       call MPI_RECV(bufftn  ,jpi_rec*jpj_rec*jpk          ,mpi_real8,idrank, 12,mpi_comm_world, status, ierr)
       call MPI_RECV(buffvatm,jpi_rec*jpj_rec              ,mpi_real8,idrank, 13,mpi_comm_world, status, ierr)
       call MPI_RECV(buffemp ,jpi_rec*jpj_rec              ,mpi_real8,idrank, 14,mpi_comm_world, status, ierr)
       call MPI_RECV(buffqsr ,jpi_rec*jpj_rec              ,mpi_real8,idrank, 15,mpi_comm_world, status, ierr)
       call MPI_RECV(buffun  ,jpi_rec*jpj_rec*jpk          ,mpi_real8,idrank, 16,mpi_comm_world, status, ierr)
       call MPI_RECV(buffvn  ,jpi_rec*jpj_rec*jpk          ,mpi_real8,idrank, 17,mpi_comm_world, status, ierr)
       call MPI_RECV(buffwn  ,jpi_rec*jpj_rec*jpk          ,mpi_real8,idrank, 18,mpi_comm_world, status, ierr)
       call MPI_RECV(buffavt ,jpi_rec*jpj_rec*jpk          ,mpi_real8,idrank, 19,mpi_comm_world, status, ierr)
       call MPI_RECV(buffe3t ,jpi_rec*jpj_rec*jpk          ,mpi_real8,idrank, 19,mpi_comm_world, status, ierr)
 ! ******* myrank 0 sets indexes of tot matrix where to place buffers of idrank
               irange    = iPe - iPd + 1
               jrange    = jPe - jPd + 1
               totistart = istart + iPd - 1
               totiend   = totistart + irange - 1
               totjstart = jstart + jPd - 1
               totjend   = totjstart + jrange - 1
               relistart = 1 + iPd - 1
               reliend   = relistart + irange - 1
               reljstart = 1 + jPd - 1
               reljend   = reljstart + jrange - 1
               do ji =totistart,totiend ! 3d vars
                    i_contribution = jpk*jpj_rec*(ji-1 - totistart+ relistart )
                     do jj =totjstart,totjend
                         j_contribution = jpk*(jj-totjstart+ reljstart-1)
                         do jk =1 , jpk
                            ind = jk + j_contribution + i_contribution
                            totsnIO (jk,jj,ji)= buffsn (ind)
                            tottnIO (jk,jj,ji)= bufftn (ind)
                            totunIO (jk,jj,ji)= buffun (ind)
                            totvnIO (jk,jj,ji)= buffvn (ind)
                            totwnIO (jk,jj,ji)= buffwn (ind)
                            totavtIO(jk,jj,ji)= buffavt(ind)
                            tote3tIO(jk,jj,ji)= buffe3t(ind)
                         enddo
                      enddo
                   enddo

             do ji =totistart,totiend  ! and 2d vars
                i_contribution = jpj_rec*(ji-1 -totistart+ relistart )
                do jj =totjstart,totjend
                   ind = (jj - totjstart+ reljstart) + i_contribution
                   totvatmIO (jj,ji)= buffvatm(ind)
                   totempIO  (jj,ji)= buffemp (ind)
                   totqsrIO  (jj,ji)= buffqsr (ind)
                enddo
             enddo
           enddo !idrank = 1, size-1
       else  ! IF LABEL 1,  if(myrank == 0)
            do ji =1 , jpi
             i_contribution= jpk*jpj * (ji - 1 )
             do jj =1 , jpj
              j_contribution=jpk*(jj-1)
              do jk =1 , jpk
                    ind         =  jk + j_contribution + i_contribution
                    buffsn (ind)= snIO (jk,jj,ji)
                    bufftn (ind)= tnIO (jk,jj,ji)
                    buffun (ind)= unIO (jk,jj,ji)
                    buffvn (ind)= vnIO (jk,jj,ji)
                    buffwn (ind)= wnIO (jk,jj,ji)
                    buffavt(ind)= avtIO(jk,jj,ji)
                    buffe3t(ind)= e3tIO(jk,jj,ji)
               enddo
              enddo
             enddo
             do ji =1 , jpi
              i_contribution=jpj * (ji-1)
              do jj =1 , jpj
                ind           = jj + i_contribution
                buffvatm (ind)= vatmIO(jj,ji)
                buffemp  (ind)= empIO (jj,ji)
                buffqsr  (ind)= qsrIO (jj,ji)
               enddo
              enddo
               call MPI_SEND(jpi  , 1,mpi_integer, 0, 1, mpi_comm_world,ierr)
               call MPI_SEND(jpj  , 1,mpi_integer, 0, 2, mpi_comm_world,ierr)
               call MPI_SEND(nimpp, 1,mpi_integer, 0, 3, mpi_comm_world,ierr)
               call MPI_SEND(njmpp, 1,mpi_integer, 0, 4, mpi_comm_world,ierr)
               call MPI_SEND(nlei , 1,mpi_integer, 0, 5, mpi_comm_world,ierr)
               call MPI_SEND(nlej , 1,mpi_integer, 0, 6, mpi_comm_world,ierr)
               call MPI_SEND(nldi , 1,mpi_integer, 0, 7, mpi_comm_world,ierr)
               call MPI_SEND(nldj , 1,mpi_integer, 0, 8, mpi_comm_world,ierr)
            call MPI_SEND(buffsn  , jpk*jpj*jpi  ,mpi_real8, 0, 11, mpi_comm_world,ierr)
            call MPI_SEND(bufftn  , jpk*jpj*jpi  ,mpi_real8, 0, 12, mpi_comm_world,ierr)
            call MPI_SEND(buffvatm, jpi*jpj      ,mpi_real8, 0, 13, mpi_comm_world,ierr)
            call MPI_SEND(buffemp , jpi*jpj      ,mpi_real8, 0, 14, mpi_comm_world,ierr)
            call MPI_SEND(buffqsr , jpi*jpj      ,mpi_real8, 0, 15, mpi_comm_world,ierr)
            call MPI_SEND(buffun  , jpk*jpj*jpi  ,mpi_real8, 0, 16, mpi_comm_world,ierr)
            call MPI_SEND(buffvn  , jpk*jpj*jpi  ,mpi_real8, 0, 17, mpi_comm_world,ierr)
            call MPI_SEND(buffwn  , jpk*jpj*jpi  ,mpi_real8, 0, 18, mpi_comm_world,ierr)
            call MPI_SEND(buffavt , jpk*jpj*jpi  ,mpi_real8, 0, 19, mpi_comm_world,ierr)
            call MPI_SEND(buffe3t , jpk*jpj*jpi  ,mpi_real8, 0, 19, mpi_comm_world,ierr)




       endif ! IF LABEL 1, if(myrank == 0)
!************* END COLLECTING DATA  *****************

! *********** START WRITING **************************

      if(myrank == 0) then ! IF LABEL 4,
         if (IsBackup) then

            call PhysDump_bkp(forcing_file, datefrom, dateTo,elapsed_time)
          else
            call PhysDump(forcing_file, datefrom, dateTo)
         endif
      endif

      endif !if ( freq_ave_phys.eq.FREQ_GROUP)



      if (WRITING_RANK_WR) tottrnIO2d = Miss_val
! ******************  DIAGNOSTIC OUTPUT   2D *******************

        
        IF (FREQ_GROUP==1)then
                elapsed_time=elapsed_time_1
                DIR='AVE_FREQ_1/'
                n_dumping_cycles=matrix_diag_2d_1_row
        end if

        if (FREQ_GROUP==2)then
                elapsed_time=elapsed_time_2
                DIR='AVE_FREQ_2/'
                n_dumping_cycles=matrix_diag_2d_2_row
        END IF


        COUNTER_VAR_2d = 1
        COUNTER_VAR_HIGH_2d = 1


        DUMPING_LOOP_2d: DO jv = 1, n_dumping_cycles

                DO ivar = 1 , nodes

                        writing_rank = writing_procs(ivar)

                        
                        IF (COUNTER_VAR_2d > JPTRA_dia_2d_wri)then
                                EXIT
                        else if (COUNTER_VAR_HIGH_2d > JPTRA_dia_2d_HIGH_wri)then
                                EXIT
                        ELSE

                                var_to_send_2D = lowfreq_table_dia_2d_wri(counter_var_2d)
                                !if (FREQ_GROUP.eq.1) var_high_to_send_2D = highfreq_table_dia_2d_wri(counter_var_high_2d)                        

        
                                if (FREQ_GROUP.eq.2) then
                                        do ji =1 , jpi
                                                i_contribution = jpj * (ji-1)
                                                do jj =1 , jpj
                                                        ind = jj + i_contribution
                                                        buffDIA2d (ind)= tra_DIA_2d_IO(var_to_send_2D,jj,ji)
                                                enddo
                                        enddo
                                else
                                        do ji =1 , jpi
                                                i_contribution = jpj * (ji-1)
                                                do jj = 1 , jpj
                                                        ind = jj + i_contribution
                                                        buffDIA2d (ind)= tra_DIA_2d_IO_high(COUNTER_VAR_HIGH_2d,jj,ji)
                                                enddo
                                        enddo
                                endif
                                counter_var_2d = counter_var_2d + 1
                                if (FREQ_GROUP.eq.1) counter_var_high_2d = counter_var_high_2d + 1

                                CALL MPI_GATHERV(buffDIA2d, sendcount_2d,MPI_DOUBLE_PRECISION, buffDIA2d_TOT,jprcv_count_2d,jpdispl_count_2d, MPI_DOUBLE_PRECISION,writing_rank, MPI_COMM_WORLD, IERR)

                        END IF
                END DO

        !---------------------------------------------------------------------------------------
        !if writitng rank assembling and dumping

                IF (WRITING_RANK_WR)then

                        ind_col = (myrank / n_ranks_per_node) +1

                        if (FREQ_GROUP.eq.2) then
                                var_to_store_diag_2d = matrix_diag_2d_2(jv,ind_col)%var_name
                        else
                                var_to_store_diag_2d = matrix_diag_2d_1(jv,ind_col)%var_name
                        end if

                        IF (var_to_store_diag_2d == "novars_input")then
                                EXIT
                        ELSE

                                do idrank = 0,mpi_glcomm_size-1
                                         irange    = iPe_a(idrank+1) - iPd_a(idrank+1) + 1
                                         jrange    = jPe_a(idrank+1) - jPd_a(idrank+1) + 1
                                         totistart = istart_a(idrank+1) + iPd_a(idrank+1) - 1
                                         totiend   = totistart + irange - 1
                                         totjstart = jstart_a(idrank+1) + jPd_a(idrank+1) - 1
                                         totjend   = totjstart + jrange - 1
                                         relistart = 1 + iPd_a(idrank+1) - 1
                                         reliend   = relistart + irange - 1
                                         reljstart = 1 + jPd_a(idrank+1) - 1
                                         reljend   = reljstart + jrange - 1
                                         do ji =totistart,totiend ! only 2d vars
                                                i_contribution = jpj_rec_a(idrank+1)*(ji-totistart+relistart -1)
                                                do jj =totjstart,totjend
                                                        ind = jj-totjstart+ reljstart +i_contribution
                                                        tottrnIO2d (jj,ji)=buffDIA2d_TOT(ind+jpdispl_count_2d(idrank+1))
                                                enddo
                                         enddo
                                enddo
                                !if (FREQ_GROUP.eq.2)write(*,*) 'CHECK ', var_to_store_diag_2d,var_to_send_2d
                                !if (FREQ_GROUP.eq.1) write(*,*)'CHECK_h', var_to_store_diag_2d, COUNTER_VAR_HIGH_2d
                                bkpname     = DIR//'ave.'//datemean//'.'//trim(var_to_store_diag_2d)//'.nc.bkp'
                                dia_file_nc = DIR//'ave.'//datemean//'.'//trim(var_to_store_diag_2d)//'.nc'

                                if (IsBackup) then
                                        CALL WRITE_AVE_2d_BKP(bkpname,var_to_store_diag_2d,datefrom, dateTo,tottrnIO2d, elapsed_time)

                                else
                                        d2f2d = REAL(tottrnIO2d(:,:),4)
                                        CALL WRITE_AVE_2d(dia_file_nc,var_to_store_diag_2d,datefrom,dateTo, d2f2d)

                                endif
                        end if
                END IF
        END DO DUMPING_LOOP_2d

        if (.not.IsBackup) then
                if (FREQ_GROUP.eq.2) then
                        tra_DIA_2d_IO(:,:,:) = 0.
                else
                        tra_DIA_2d_IO_HIGH(:,:,:) = 0.
                endif
        endif

!------------------------------------------------------------------------------------------------

        if (WRITING_RANK_WR) tottrnIO = Miss_val
! ! ******************  3D DIAGNOSTIC OUTPUT   *******************

        IF (FREQ_GROUP==1)then
                elapsed_time=elapsed_time_1
                DIR='AVE_FREQ_1/'
                !matrix_col=nodes
                n_dumping_cycles=matrix_diag_1_row
        end if

        if (FREQ_GROUP==2)then
                elapsed_time=elapsed_time_2
                DIR='AVE_FREQ_2/'
                !matrix_col=nodes
                n_dumping_cycles=matrix_diag_2_row
        END IF


        COUNTER_VAR_diag = 1
        COUNTER_VAR_diag_HIGH = 1


        DUMPING_LOOP_3d: DO jv = 1, n_dumping_cycles

                DO ivar = 1 , nodes

                        writing_rank = writing_procs(ivar)


                        IF (COUNTER_VAR_diag > JPTRA_dia_wri)then
                                EXIT
                        else if (COUNTER_VAR_diag_HIGH > JPTRA_dia_HIGH_wri)then
                                EXIT
                        ELSE

                                var_to_send = lowfreq_table_dia_wri(counter_var_diag)
                                !if (FREQ_GROUP.eq.1) var_high_to_send = highfreq_table_dia_wri(counter_var_diag_high)

                                if (FREQ_GROUP.eq.2) then
                                        do ji =1, jpi
                                                i_contribution= jpk*jpj * (ji - 1 )
                                                do jj =1 , jpj
                                                        j_contribution=jpk*(jj-1)
                                                        do jk =1 , jpk
                                                                ind = jk + j_contribution + i_contribution
                                                                buffDIA(ind) = tra_DIA_IO(var_to_send, jk,jj,ji)
                                                        enddo
                                                enddo
                                        enddo
                                else
                                        do ji =1, jpi
                                                i_contribution= jpk*jpj * (ji - 1 )
                                                do jj =1 , jpj
                                                        j_contribution=jpk*(jj-1)
                                                        do jk =1 , jpk
                                                                ind = jk + j_contribution + i_contribution
                                                                buffDIA(ind) = tra_DIA_IO_HIGH(COUNTER_VAR_diag_HIGH, jk,jj,ji)
                                                        enddo
                                                enddo
                                        enddo
                                end if
                                !if (FREQ_GROUP.eq.1)write(*,*)'CHECK_h_before', COUNTER_VAR_diag_HIGH
                                counter_var_diag = counter_var_diag + 1
                                if (FREQ_GROUP.eq.1) counter_var_diag_high = counter_var_diag_high + 1

                                !GATHERV TO THE WRITING RANK

                                CALL MPI_GATHERV(buffDIA, sendcount,MPI_DOUBLE_PRECISION, buffDIA_TOT,jprcv_count, jpdispl_count,MPI_DOUBLE_PRECISION, writing_rank,MPI_COMM_WORLD, IERR)
                        END IF
                END DO

! *********** START WRITING **************************
                IF (WRITING_RANK_WR)then

                        ind_col = (myrank / n_ranks_per_node)+1

                        if (FREQ_GROUP.eq.2) then
                                var_to_store_diag = matrix_diag_2(jv,ind_col)%var_name
                        else
                                var_to_store_diag = matrix_diag_1(jv,ind_col)%var_name
                        end if

                        IF (var_to_store_diag == "novars_input")then
                                EXIT
                        ELSE

                                do idrank = 0,mpi_glcomm_size-1
                                        irange    = iPe_a(idrank+1) - iPd_a(idrank+1) + 1
                                        jrange    = jPe_a(idrank+1) - jPd_a(idrank+1) + 1
                                        totistart = istart_a(idrank+1) + iPd_a(idrank+1) - 1
                                        totiend   = totistart + irange - 1
                                        totjstart = jstart_a(idrank+1) + jPd_a(idrank+1) - 1
                                        totjend   = totjstart + jrange - 1
                                        relistart = 1 + iPd_a(idrank+1) - 1
                                        reliend   = relistart + irange - 1
                                        reljstart = 1 + jPd_a(idrank+1) - 1
                                        reljend   = reljstart + jrange - 1

                                        do ji =totistart,totiend ! 3d vars
                                                i_contribution = jpk*jpj_rec_a(idrank+1)*(ji-1 -totistart+relistart )
                                                do jj =totjstart,totjend
                                                        j_contribution = jpk*(jj-totjstart+ reljstart-1)
                                                        do jk =1 , jpk
                                                                ind = jk + j_contribution + i_contribution
                                                                tottrnIO(jk,jj,ji)= buffDIA_TOT (ind+jpdispl_count(idrank+1))
                                                        enddo
                                                enddo
                                        enddo
                                enddo
                                !if (FREQ_GROUP.eq.2)write(*,*) 'CHECK ', var_to_store_diag, var_to_send
                                !if (FREQ_GROUP.eq.1)write(*,*) 'CHECK_h', var_to_store_diag, COUNTER_VAR_diag_HIGH
                                bkpname     = DIR//'ave.'//datemean//'.'//trim(var_to_store_diag)//'.nc.bkp'
                                dia_file_nc = DIR//'ave.'//datemean//'.'//trim(var_to_store_diag)//'.nc'
              
                                if (IsBackup) then
                                        CALL WRITE_AVE_BKP(bkpname,var_to_store_diag,datefrom, dateTo,tottrnIO,elapsed_time,deflate_ave, deflate_level_ave)
                                else
                                        CALL WRITE_AVE(dia_file_nc,var_to_store_diag,datefrom,dateTo, tottrnIO,deflate_ave, deflate_level_ave)
                                endif


                        END IF
                END IF
        END DO DUMPING_LOOP_3d

        if (.not.IsBackup) then
                if (FREQ_GROUP.eq.2) then
                        tra_DIA_IO(:,:,:,:) = 0.
                else
                        tra_DIA_IO_HIGH(:,:,:,:) = 0.
                endif
        endif

        if ((.not.IsBackup).and.( freq_ave_phys.eq.FREQ_GROUP) ) then

                snIO     = 0.
                tnIO     = 0.
                vatmIO   = 0.
                empIO    = 0.
                qsrIO    = 0.
                unIO     = 0.
                vnIO     = 0.
                wnIO     = 0.
                avtIO    = 0.
                e3tIO    = 0
        endif


        end SUBROUTINE diadump
