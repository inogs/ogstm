      SUBROUTINE read_float(filename, var_to_store, shift, reading_proc,matrix_scatter_local)
      USE calendar
      USE myalloc
      USE TIME_MANAGER
      USE ogstm_mpi_module

      USE MPI_GATHER_INFO
      USE MATRIX_VARS
      USE NODES_MODULE
      USE DTYPE_PROCS_STRING_MODULE

      !character(LEN=30) nomefile
      character, INTENT(IN) :: filename*(*) ,var_to_store*(*)
      INTEGER, INTENT(IN) :: shift, reading_proc
      double precision,INTENT(OUT) :: matrix_scatter_local(jpk,jpj,jpi)

      double precision, dimension(jpk,jpjglo,jpiglo) :: buf_glo

      INTEGER idrank, ierr, istart, jstart, iPe, iPd, jPe, jPd, status(MPI_STATUS_SIZE)
      INTEGER irange, jrange
      INTEGER totistart, totiend, relistart, reliend
      INTEGER totjstart, totjend, reljstart, reljend
      INTEGER jk,jj,ji,i,j,k
      INTEGER ind,ind_loc, i_contribution, j_contribution

      DOUBLE PRECISION, allocatable, dimension(:) :: buff_scatterv_tot
      double precision, allocatable :: buff_scatter(:)
      
      CHARACTER(LEN=20)  DIR,date_from,date_To
      CHARACTER(LEN=56) output_file_nc

      allocate(buff_scatter (jpi* jpj* jpk))


      buff_scatter = huge(buff_scatter(1))
      
      ALLOCATE (buff_scatterv_tot(total_dim_jpijk))
      buff_scatterv_tot = huge(buff_scatterv_tot(1))

      call mppsync()

      if(myrank.eq.0) then
              call readnc_slice_float_INITIAL(filename,var_to_store,buf_glo,ingv_lon_shift)

      call mppsync()

      !!! il proc 0 deve trasformare la matrice 3d in un array monodimensionale, secondo gli indici dei processori

      if(lwp) then
        ind=1
        DO idrank = 0,mpi_glcomm_size-1
                totistart = istart_a(idrank + 1) 
                totiend   = istart_a(idrank + 1) + jpi_rec_a(idrank + 1) -1
                totjstart = jstart_a(idrank + 1)
                totjend   = jstart_a(idrank + 1) + jpj_rec_a(idrank + 1) -1
        
                do ji =totistart,totiend
                        do jj =totjstart,totjend
                                do jk =1, jpk
                                        buff_scatterv_tot(ind)=buf_glo(jk,jj,ji)
                                        ind=ind+1
                                enddo
                        enddo
                enddo
        END DO
      end if
      

      !adesso il proc 0 manda i pezzi del buff scatterv tot a tutti i procs
      if(myrank.eq.0) write(*,*) 'assembling done'

      call mppsync()

      CALL MPI_SCATTERV(buff_scatterv_tot, jprcv_count, jpdispl_count, MPI_DOUBLE_PRECISION, buff_scatter, sendcount,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERR)
      
      if(myrank.eq.0) write(*,*) 'scatterv done'

      !iogni processore riceve il buffer e si costruisce il suo pezzo di matrice

      call mppsync()

      ind_loc = 1
      do ji =1 , jpi
          do jj =1 , jpj
                  do jk =1 , jpk
                          matrix_scatter_local(jk,jj,ji)=buff_scatter(ind_loc)
                          ind_loc = ind_loc + 1
                  enddo
          enddo
      enddo
      
      if(lwp) write(*,*) "reconstruction done"


      END SUBROUTINE read_float
