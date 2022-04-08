      SUBROUTINE read_float(filename, var_to_store, buf, shift, reading_proc)
      USE calendar
      USE myalloc
      USE TIME_MANAGER
      USE ogstm_mpi_module

      USE MPI_GATHER_INFO
      USE MATRIX_VARS
      USE NODES_MODULE
      USE DTYPE_PROCS_STRING_MODULE

      character(LEN=30) nomefile
      character, INTENT(IN) :: filename*(*) ,var_to_store*(*)

      CHARACTER(LEN=20)  var_to_store
      double precision, dimension(jpk,jpjglo,jpiglo) :: buf_glo, test_buf_glo

      INTEGER idrank, ierr, istart, jstart, iPe, iPd, jPe, jPd, status(MPI_STATUS_SIZE)
      INTEGER irange, jrange
      INTEGER totistart, totiend, relistart, reliend
      INTEGER totjstart, totjend, reljstart, reljend
      INTEGER jk,jj,ji,i,j,k
      INTEGER ind, i_contribution, j_contribution

      DOUBLE PRECISION, allocatable, dimension(:) :: buff_scatterv_tot, test_buff_scatterv_tot
      integer, dimension(mpi_glcomm_size)::sendcount_array
      double precision, allocatable :: buff_scatter(:), test_buff_scatter(:)
      double precision, dimension(jpk,jpj,jpi) :: matrix_scatter_local

      allocate(buff_scatter (jpi_max* jpj_max* jpk))
      allocate(test_buff_scatter (jpi_max* jpj_max* jpk))


      buff_scatter = huge(buff_scatter(1))
      test_buff_scatter = huge(test_buff_scatter(1))
      
      ALLOCATE (buff_scatterv_tot(jpiglo* jpjglo* jpk))
              ALLOCATE (test_buff_scatterv_tot(jpiglo* jpjglo* jpk))
              !ALLOCATE (sendcount_array(mpi_glcomm_size))
              buff_scatterv_tot = huge(buff_scatterv_tot(1))
              test_buff_scatterv_tot = huge(test_buff_scatterv_tot(1))
              if(lwp) write(*,*) 'allocation done'


      call mppsync()
      CALL MPI_Gather(sendcount, 1, MPI_INT,sendcount_array, 1, MPI_INT, 0,MPI_COMM_WORLD, IERR)
      write(*,*) "gather done"

      if (lwp) then
        DO idrank = 0,mpi_glcomm_size-1
                write(*,*) 'rank ', idrank+1, 'sendcount is ', sendcount_array(idrank+1)
        END DO
      end if


      call mppsync()



      nomefile='FORCINGS/U19951206-12:00:00.nc'

! Starting I/O
! U  *********************************************************
      nomefile = 'FORCINGS/U'//datestring//'.nc'
      if(lwp) write(*,'(A,I4,A,A)') "LOAD_PHYS --> I am ", myrank, " starting reading forcing fields from ", nomefile(1:30)
      if(lwp) then
              call readnc_slice_float_INITIAL(nomefile,'vozocrtx',buf_glo,ingv_lon_shift)
              udta(:,:,:,2) = buf * umask
              DIR='AVE_FREQ_1/'
              var_to_store='vozocrtx'
              output_file_nc=trim(DIR)//'ave.TEST.'//trim(var_to_store)//'.nc'
              date_from='20220302-12:00:00'
              date_To='20220302-13:00:00'


      CALL WRITE_AVE(output_file_nc,var_to_store,date_from, date_To, buf_glo, deflate_ave, deflate_level_ave)
      end if

      call mppsync()
      if (lwp) write(*,*) 'write ave done'
     ! STOP

      !!! il proc 0 deve trasformare la matrice 3d in un array monodimensionale, secondo gli indici dei processori

      if(lwp) then

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
                                         buff_scatterv_tot(ind+jpdispl_count(idrank+1))=buf_glo(jk,jj,ji)
                                enddo
                        enddo
                enddo
        END DO

     end if

    !sendcount_array
     !mpi gather
     !ttui devono mandare il loro sendcount a 0
     !   call mpi_scatterv()

        !adesso il proc 0 manda i pezzi del buff scatterv tot a tutti i procs
     if(lwp) write(*,*) 'assembling done'

     call mppsync()

     CALL MPI_SCATTERV(buff_scatterv_tot,sendcount_array, jpdispl_count, MPI_DOUBLE_PRECISION, buff_scatter, sendcount,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERR)


     call mppsync()
     if(lwp) write(*,*) "scatterv done"

     !iogni processore rieve il buffer e si costruisce il suo pezzo di matrice

     do ji =1 , jpi
         i_contribution= jpk*jpj * (ji - 1 )
         do jj =1 , jpj
                 j_contribution=jpk*(jj-1)
                 do jk =1 , jpk
                         ind =  jk + j_contribution + i_contribution
                         matrix_scatter_local(jk,jj,ji)=buff_scatter(ind)
                 enddo
         enddo
     enddo

     call mppsync()
     if(lwp) write(*,*) "reconstruction done"












