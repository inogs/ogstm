SUBROUTINE trcdit(datemean,datefrom,dateTo,FREQ_GROUP)
        !---------------------------------------------------------------------
        !
        !                       ROUTINE trcdit
        !
        !                     ******************
        !  gcoidess develop
        !
        !  Purpose :
        !  ---------
        !     Standard output of passive tracer : concentration fields



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

        ! local declarations
        ! ==================
        CHARACTER(LEN=17), INTENT(IN) :: datemean, dateFrom, dateTo
        INTEGER, INTENT(IN) :: FREQ_GROUP ! 1 = HIGH FREQ, 2 = LOW FREQ #cannot change value inside te value
        double precision ::  Miss_val =1.e20

        INTEGER jk,jj,ji,i,j,k
        INTEGER ind, i_contribution, j_contribution
        double precision :: elapsed_time

        !new declarations
        INTEGER counter_var, counter_var_high, new_counter_var, new_counter_var_high, nVARS, jv, ivar, n_dumping_cycles
        INTEGER col_var, row_var, writing_rank
        CHARACTER(len=20), DIMENSION(nodes) :: matrix_row_to_write

        CHARACTER(LEN=56) output_file_nc  ! AVE_FREQ_1/ave.20091231-12:00:00.P1n.nc
        CHARACTER(LEN=20) var
        CHARACTER(LEN=60) bkpname
        CHARACTER(LEN=11) DIR
        logical IsBackup

        CHARACTER(LEN=20)  var_to_store

        INTEGER :: ind_col
        DOUBLE PRECISION :: start_time_trcdit_info,finish_time_trcdit_info, proctime_time_trcdit_info, max_time_trcdit_info
        DOUBLE PRECISION :: gatherv_fin_time,gatherv_init_time,gatherv_delta_time,gatherv_sum_time,gatherv_mean_time
        DOUBLE PRECISION :: writing_rank_fin_time, writing_rank_init_time,writing_rank_delta_time, writing_rank_sum_time


        !----------------------------------------------------------------------
        ! statement functions
        ! ===================


        INTEGER idrank, ierr, istart, jstart, iPe, iPd, jPe, jPd, status(MPI_STATUS_SIZE)
        INTEGER irange, jrange
        INTEGER totistart, totiend, relistart, reliend
        INTEGER totjstart, totjend, reljstart, reljend

        ! ----------------------------------------
        IsBackup =  (datemean.eq.dateTo)
        if (lwp) write(*,*) 'trcdit IsBackup = ',IsBackup, ' group ' ,FREQ_GROUP
        ! ----------------------------------------
        bkpname  = 'AVE_FREQ_1/ave.20111231-15:30:00.N1p.nc.bkp'
        call mppsync()


        SELECT CASE (FREQ_GROUP)
        CASE (1) 
                elapsed_time=elapsed_time_1
                DIR='AVE_FREQ_1/'
                n_dumping_cycles = matrix_state_1_row
        CASE (2) 
                elapsed_time=elapsed_time_2
                DIR='AVE_FREQ_2/'
                n_dumping_cycles = matrix_state_2_row
        END SELECT

        if (WRITING_RANK_WR)then
                tottrnIO = Miss_val
        endif
      
        !-----------------------
        !starting big loop
        start_time_trcdit_info= MPI_Wtime()

        COUNTER_VAR = 1
        COUNTER_VAR_HIGH = 1

        DUMPING_LOOP: DO jv = 1, n_dumping_cycles

                gatherv_init_time = MPI_Wtime()

                DO ivar = 1 , nodes!number of variables for each round corresponds to the number of nodes

                        writing_rank = writing_procs(ivar)

                        IF (COUNTER_VAR > JPTRA)then
                                EXIT
                        else if (COUNTER_VAR_HIGH > JPTRA_HIGH)then
                                EXIT
                        ELSE

                                if (FREQ_GROUP.eq.2) then
                                        do ji =1 , jpi
                                                i_contribution= jpk*jpj * (ji - 1 )
                                                do jj =1 , jpj
                                                        j_contribution=jpk*(jj-1)
                                                        do jk =1 , jpk
                                                                ind =  jk + j_contribution + i_contribution
                                                                bufftrn   (ind)= traIO( jk,jj,ji,counter_var)
                                                        enddo
                                                enddo
                                        enddo
                                else ! FREQ_GROUP.eq.1
                                        do ji =1 , jpi
                                                i_contribution= jpk*jpj * (ji - 1 )
                                                do jj =1 , jpj
                                                        j_contribution=jpk*(jj-1)
                                                        do jk =1 , jpk
                                                                ind =  jk + j_contribution + i_contribution
                                                                bufftrn   (ind)= traIO_HIGH(jk,jj,ji,counter_var_high)
                                                        enddo
                                                enddo
                                        enddo
                                endif

                                counter_var = counter_var + 1
                                if (FREQ_GROUP.eq.1) counter_var_high = counter_var_high + 1

                                !GATHERV TO THE WRITING RANK
                                CALL MPI_GATHERV(bufftrn, sendcount, MPI_DOUBLE_PRECISION, bufftrn_TOT, jprcv_count,& 
                                                jpdispl_count, MPI_DOUBLE_PRECISION, writing_rank, MPI_COMM_WORLD, IERR)

                        END IF

                END DO
                gatherv_fin_time = MPI_Wtime()

                gatherv_delta_time = gatherv_fin_time - gatherv_init_time
                CALL MPI_Reduce( gatherv_delta_time, gatherv_sum_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD,IERROR)


                ! *************** COLLECTING DATA *****************************

                IF (WRITING_RANK_WR)then

                        writing_rank_init_time = MPI_Wtime()

                        ind_col = (myrank / n_ranks_per_node)+1

                        if (FREQ_GROUP.eq.2) then
                                var_to_store = matrix_state_2(jv,ind_col)%var_name
                        else
                                var_to_store = matrix_state_1(jv,ind_col)%var_name
                        end if

                        IF (var_to_store == "novars_input")then
                                EXIT
                        ELSE
                                DO idrank = 0,mpi_glcomm_size-1

                                        ! ******* WRITING RANK sets indexes of tot matrix where to place buffers of idrank
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

                                        ! **** ASSEMBLING *** WRITING RANK  puts in tot matrix buffer received by idrank
                                        do ji =totistart,totiend
                                         i_contribution   = jpk*jpj_rec_a(idrank+1)*(ji-1-totistart+ relistart)
                                          do jj =totjstart,totjend
                                           j_contribution = jpk*(jj-1-totjstart+ reljstart)
                                           do jk =1, jpk
                                            ind = jk + j_contribution + i_contribution
                                            tottrnIO(jk,jj,ji)= bufftrn_TOT(ind+jpdispl_count(idrank+1))
                                           enddo
                                          enddo
                                        enddo
                                END DO

                                output_file_nc = DIR//'ave.'//datemean//'.'//trim(var_to_store)//'.nc'
                                bkpname = DIR//'ave.'//datemean//'.'//trim(var_to_store)//'.nc.bkp'

                                if (IsBackup) then
                                        CALL WRITE_AVE_BKP(bkpname,var_to_store,datefrom,&
                                                dateTo,tottrnIO,elapsed_time, deflate_ave, deflate_level_ave)
                                else
                                        CALL WRITE_AVE(output_file_nc,var_to_store,datefrom,&
                                                dateTo, tottrnIO, deflate_ave, deflate_level_ave)
                                endif
                        END IF
                        !writing_rank_fin_time = MPI_Wtime()
                        !writing_rank_delta_time = writing_rank_fin_time - writing_rank_init_time
                        !writing_rank_sum_time = writing_rank_delta_time + writing_rank_sum_time
                        !write(*,*)'writingtime', writing_rank_sum_time,'   ',jv, '   ', myrank
                END IF
        END DO DUMPING_LOOP

        if (.not.IsBackup) then
                if (FREQ_GROUP.eq.2) then
                        traIO(:,:,:,:) = 0.  !      we reset matrix for new average
                else
                        traIO_HIGH(:,:,:,:) = 0.
                endif
        end if

        finish_time_trcdit_info= MPI_Wtime()
        proctime_time_trcdit_info=finish_time_trcdit_info - start_time_trcdit_info

        CALL MPI_Reduce( proctime_time_trcdit_info,max_time_trcdit_info, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD,IERROR)

        !if(myrank == 0) then
        !        write(*,*) 'TRCDIT TIME is', max_time_trcdit_info
        !end if


END SUBROUTINE trcdit
