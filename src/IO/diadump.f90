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
      CHARACTER(LEN=56) phys_file_nc
      CHARACTER(LEN=20)  var

      INTEGER idrank, ierr, istart, jstart, iPe, iPd, jPe, jPd, status(MPI_STATUS_SIZE)
      INTEGER irange, jrange
      INTEGER totistart, totiend, relistart, reliend
      INTEGER totjstart, totjend, reljstart, reljend
      double precision ::  Miss_val =1.e20
      INTEGER :: nVars, counter_var_2d, counter_var_high_2d,counter_var_diag, counter_var_diag_high
      INTEGER :: counter_var_phys_2d,counter_var_phys_high_2d,counter_var_phys, counter_var_phys_high
      CHARACTER(LEN=20) ::  var_to_store_diag_2d, var_to_store_diag
      CHARACTER(LEN=20) ::  var_to_store_phys_2d, var_to_store_phys
      INTEGER :: n_dumping_cycles, jv, ivar, writing_rank, ind_col
      INTEGER :: var_to_send_2D, var_high_to_send_2D
      INTEGER :: var_to_send, var_high_to_send
     ! call mppsync()
! ----------------------------------------
      IsBackup =  (datemean.eq.dateTo)
      if (lwp) write(*,*) 'diadump IsBackup = ',IsBackup, ' group ' ,FREQ_GROUP
! ----------------------------------------

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

                                CALL MPI_GATHERV(buffDIA2d,sendcount_2d,MPI_DOUBLE_PRECISION, &
                                                 buffDIA2d_TOT,jprcv_count_2d,jpdispl_count_2d,&
                                                 MPI_DOUBLE_PRECISION,writing_rank, MPI_COMM_WORLD, IERR)

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
                                        CALL WRITE_AVE_2d_BKP(bkpname,var_to_store_diag_2d,datefrom, &
                                                              dateTo,tottrnIO2d, elapsed_time)

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

                                CALL MPI_GATHERV(buffDIA, sendcount,MPI_DOUBLE_PRECISION,buffDIA_TOT,jprcv_count, &
                                    jpdispl_count,MPI_DOUBLE_PRECISION, writing_rank,MPI_COMM_WORLD, IERR)
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
                                        CALL WRITE_AVE_BKP(bkpname,var_to_store_diag,datefrom, dateTo,&
                                                       tottrnIO,elapsed_time,deflate_ave, deflate_level_ave)
                                else
                                        CALL WRITE_AVE(dia_file_nc,var_to_store_diag,datefrom,dateTo, &
                                         tottrnIO,deflate_ave, deflate_level_ave)
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
!-------------------------------------------------
        ! ****************** PHYSC OUTPUT   2D *******************
        




        if(freq_ave_phys.eq.FREQ_GROUP ) then

        tra_PHYS_2d_IO(1,:,:) = vatmIO
        tra_PHYS_2d_IO(2,:,:) = empIO
        tra_PHYS_2d_IO(3,:,:) = qsrIO

        tra_PHYS_2d_IO_high(1,:,:) = vatmIO
        tra_PHYS_2d_IO_high(2,:,:) = empIO
        !write(*,*) 'copy is',tra_PHYS_2d_IO_high(2,5,5)
        tra_PHYS_2d_IO_high(3,:,:) = qsrIO
        !write(*,*) 'copy is',tra_PHYS_2d_IO_high(3,5,5)


        IF (freq_ave_phys==1)then
                elapsed_time=elapsed_time_1
                DIR='AVE_FREQ_1/'
                n_dumping_cycles=matrix_phys_2d_1_row
        end if

        if (freq_ave_phys==2)then
                elapsed_time=elapsed_time_2
                DIR='AVE_FREQ_2/'
                n_dumping_cycles=matrix_phys_2d_2_row
        END IF


        COUNTER_VAR_phys_2d = 1
        COUNTER_VAR_phys_HIGH_2d = 1


        DUMPING_LOOP_2d_phys: DO jv = 1, n_dumping_cycles

                DO ivar = 1 , nodes

                        writing_rank = writing_procs(ivar)

                        !write(*,*)'phys 2d wri number' ,JPTRA_phys_2d_HIGH_wri
                        IF (freq_ave_phys==2 .and. COUNTER_VAR_phys_2d > JPTRA_phys_2d_wri)then
                                EXIT
                        else if (freq_ave_phys==1 .and. COUNTER_VAR_phys_HIGH_2d > JPTRA_phys_2d_HIGH_wri)then
                                EXIT
                        ELSE
                                if(freq_ave_phys==2) then
                                        var_to_send_2D = lowfreq_table_phys_2d_wri(counter_var_phys_2d)
                                else
                                        var_to_send_2D = highfreq_table_phys_2d_wri(counter_var_phys_high_2d)
                                        !write(*,*) 'var to send 2d is', var_to_send_2D
                                end if

                                if (freq_ave_phys.eq.2) then
                                        do ji =1 , jpi
                                                i_contribution = jpj * (ji-1)
                                                do jj =1 , jpj
                                                        ind = jj +i_contribution
                                                        buffPHYS2d (ind)=tra_PHYS_2d_IO(var_to_send_2D,jj,ji)
                                                enddo
                                        enddo
                                else
                                        do ji =1 , jpi
                                                i_contribution = jpj * (ji-1)
                                                do jj = 1 , jpj
                                                        ind = jj +i_contribution
                                                        buffPHYS2d (ind)=tra_PHYS_2d_IO_high(var_to_send_2D,jj,ji)
                                                enddo
                                        enddo
                                        !write(*,*) 'valeu in buffer is', tra_PHYS_2d_IO_high(var_to_send_2D,5,5)
                                        !write(*,*) &
                                     !'valeu in buffer is,second print',tra_PHYS_2d_IO_high(3,5,5)
                                endif
                                counter_var_phys_2d = counter_var_phys_2d + 1
                                if (freq_ave_phys.eq.1) counter_var_phys_high_2d = counter_var_phys_high_2d + 1

                                CALL MPI_GATHERV(buffPHYS2d,sendcount_2d,MPI_DOUBLE_PRECISION,buffPHYS2d_TOT,jprcv_count_2d, &
                                jpdispl_count_2d,MPI_DOUBLE_PRECISION,writing_rank, MPI_COMM_WORLD, IERR)

                        END IF
                END DO

        !---------------------------------------------------------------------------------------
        !if writitng rank assembling and dumping

                IF (WRITING_RANK_WR)then

                        ind_col = (myrank / n_ranks_per_node) +1

                        if (freq_ave_phys.eq.2) then
                                var_to_store_phys_2d = matrix_phys_2d_2(jv,ind_col)%var_name
                        else
                                var_to_store_phys_2d = matrix_phys_2d_1(jv,ind_col)%var_name
                        end if

                        IF (var_to_store_phys_2d == "novars_input")then
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
                                                        tottrnIO2d(jj,ji)=buffPHYS2d_TOT(ind+jpdispl_count_2d(idrank+1))
                                                enddo
                                         enddo
                                enddo
                                bkpname     =DIR//'ave.'//datemean//'.'//trim(var_to_store_phys_2d)//'.nc.bkp'
                                phys_file_nc =DIR//'ave.'//datemean//'.'//trim(var_to_store_phys_2d)//'.nc'

                                if (IsBackup) then
                                        CALL WRITE_AVE_2d_BKP(bkpname,var_to_store_phys_2d,datefrom, &
                                        dateTo,tottrnIO2d, elapsed_time)

                                else
                                        d2f2d = REAL(tottrnIO2d(:,:),4)
                                        CALL WRITE_AVE_2d(phys_file_nc,var_to_store_phys_2d,datefrom,&
                                                         dateTo, d2f2d)

                                endif
                        end if
                END IF
        END DO DUMPING_LOOP_2d_phys

        if (.not.IsBackup) then
                if (freq_ave_phys.eq.2) then
                        tra_PHYS_2d_IO(:,:,:) = 0.
                else
                        tra_PHYS_2d_IO_HIGH(:,:,:) = 0.
                endif
        endif


        if (WRITING_RANK_WR) tottrnIO = Miss_val
!-------------------------------------------------------------------------


       ! ! ******************  3D PHYS OUTPUT   *******************
        
        tra_PHYS_IO(1,:,:,:) = snIO
        tra_PHYS_IO(2,:,:,:) = tnIO
        tra_PHYS_IO(3,:,:,:) = wnIO
        tra_PHYS_IO(4,:,:,:) = avtIO
        tra_PHYS_IO(5,:,:,:) = e3tIO
        tra_PHYS_IO(6,:,:,:) = unIO
        tra_PHYS_IO(7,:,:,:) = vnIO

        tra_PHYS_IO_high(1,:,:,:) = snIO
        tra_PHYS_IO_high(2,:,:,:) = tnIO
        tra_PHYS_IO_high(3,:,:,:) = wnIO
        tra_PHYS_IO_high(4,:,:,:) = avtIO
        tra_PHYS_IO_high(5,:,:,:) = e3tIO
        tra_PHYS_IO_high(6,:,:,:) = unIO
        tra_PHYS_IO_high(7,:,:,:) = vnIO

        IF (freq_ave_phys==1)then
                elapsed_time=elapsed_time_1
                DIR='AVE_FREQ_1/'
                !matrix_col=nodes
                n_dumping_cycles=matrix_phys_1_row
        end if

        if (freq_ave_phys==2)then
                elapsed_time=elapsed_time_2
                DIR='AVE_FREQ_2/'
                !matrix_col=nodes
                n_dumping_cycles=matrix_phys_2_row
        END IF


        COUNTER_VAR_phys = 1
        COUNTER_VAR_phys_HIGH = 1


        DUMPING_LOOP_3d_phys: DO jv = 1, n_dumping_cycles

                DO ivar = 1 , nodes

                        writing_rank = writing_procs(ivar)

                        IF (freq_ave_phys==2 .and. COUNTER_VAR_phys > JPTRA_phys_wri)then
                                EXIT
                        else if (freq_ave_phys==1 .and. COUNTER_VAR_phys_HIGH > JPTRA_phys_HIGH_wri)then
                                EXIT
                        ELSE
                       
                                if(freq_ave_phys==2) then
                                        var_to_send = lowfreq_table_phys_wri(counter_var_phys)
                                else
                                        var_to_send = highfreq_table_phys_wri(counter_var_phys_high)
                                end if

                                if (freq_ave_phys.eq.2) then
                                        do ji =1, jpi
                                                i_contribution= jpk*jpj * (ji -1 )
                                                do jj =1 , jpj
                                                        j_contribution=jpk*(jj-1)
                                                        do jk =1 , jpk
                                                                ind = jk + j_contribution + i_contribution
                                                                buffPHYS(ind) = tra_PHYS_IO(var_to_send, jk,jj,ji)
                                                        enddo
                                                enddo
                                        enddo
                                else
                                        do ji =1, jpi
                                                i_contribution= jpk*jpj * (ji -1 )
                                                do jj =1 , jpj
                                                        j_contribution=jpk*(jj-1)
                                                        do jk =1 , jpk
                                                                ind = jk + j_contribution + i_contribution
                                                                buffPHYS(ind) = tra_PHYS_IO_HIGH(var_to_send, jk,jj,ji)
                                                        enddo
                                                enddo
                                        enddo
                                end if
                                counter_var_phys = counter_var_phys + 1
                                if (freq_ave_phys.eq.1) counter_var_phys_high =counter_var_phys_high + 1

                                !GATHERV TO THE WRITING RANK

                                CALL MPI_GATHERV(buffPHYS,sendcount,MPI_DOUBLE_PRECISION,buffPHYS_TOT,jprcv_count, &
                                                 jpdispl_count,MPI_DOUBLE_PRECISION, writing_rank,MPI_COMM_WORLD, IERR)
                        END IF
                END DO

! *********** START WRITING **************************
                IF (WRITING_RANK_WR)then

                        ind_col = (myrank / n_ranks_per_node)+1

                        if (freq_ave_phys.eq.2) then
                                var_to_store_phys = matrix_phys_2(jv,ind_col)%var_name
                        else
                                var_to_store_phys = matrix_phys_1(jv,ind_col)%var_name
                        end if

                        IF (var_to_store_phys == "novars_input")then
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
                                                                tottrnIO(jk,jj,ji)= buffPHYS_TOT (ind+jpdispl_count(idrank+1))
                                                        enddo
                                                enddo
                                        enddo
                                enddo
                                bkpname     =DIR//'ave.'//datemean//'.'//trim(var_to_store_phys)//'.nc.bkp'
                                phys_file_nc =DIR//'ave.'//datemean//'.'//trim(var_to_store_phys)//'.nc'

                                if (IsBackup) then
                                        CALL WRITE_AVE_BKP(bkpname,var_to_store_phys,datefrom,dateTo,tottrnIO,& 
                                             elapsed_time,deflate_ave, deflate_level_ave)
                                else
                                        CALL WRITE_AVE(phys_file_nc,var_to_store_phys,datefrom,dateTo, tottrnIO, & 
                                                   deflate_ave,deflate_level_ave)
                                endif


                        END IF
                END IF
        END DO DUMPING_LOOP_3d_phys

        if (.not.IsBackup) then
                if (freq_ave_phys.eq.2) then
                        tra_PHYS_IO(:,:,:,:) = 0.
                else
                        tra_PHYS_IO_HIGH(:,:,:,:) = 0.
                endif
        endif

        end if
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
