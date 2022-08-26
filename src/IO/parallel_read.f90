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
      real, dimension(jpi,jpj,jpk) :: copy_in

      real, dimension(jpiglo,jpjglo,jpk) :: buf_glo

      INTEGER idrank, ierr, istart, jstart, status(MPI_STATUS_SIZE)
      INTEGER totistart, totiend 
      INTEGER totjstart, totjend
      INTEGER jk,jj,ji,i,j,k
      INTEGER ind,ind_loc

      real, allocatable, dimension(:) :: buff_scatterv_tot
      real, allocatable :: buff_scatter(:)

      allocate(buff_scatter (jpi* jpj* jpk))
      buff_scatter = huge(buff_scatter(1))
      ALLOCATE (buff_scatterv_tot(total_dim_jpijk))
      buff_scatterv_tot = huge(buff_scatterv_tot(1))

      if(myrank.eq.0) then
              call readnc_float_global(filename,var_to_store,buf_glo,ingv_lon_shift)
      end if

      !!! il proc 0 deve trasformare la matrice 3d in un array monodimensionale, secondo gli indici dei processori

      if(myrank.eq.0) then
        ind=1
        DO idrank = 0,mpi_glcomm_size-1
                totistart = istart_a(idrank + 1) 
                totiend   = istart_a(idrank + 1) + jpi_rec_a(idrank + 1) -1
                totjstart = jstart_a(idrank + 1)
                totjend   = jstart_a(idrank + 1) + jpj_rec_a(idrank + 1) -1
        
                do jk =1, jpk
                        do jj =totjstart,totjend
                                do ji =totistart,totiend
                                        buff_scatterv_tot(ind)=buf_glo(ji,jj,jk)
                                        ind=ind+1
                                enddo
                        enddo
                enddo
        END DO
      end if
      

      !adesso il proc 0 manda i pezzi del buff scatterv tot a tutti i procs

      call mppsync()

      CALL MPI_SCATTERV(buff_scatterv_tot, jprcv_count, jpdispl_count, MPI_REAL, buff_scatter, sendcount,MPI_REAL, 0, MPI_COMM_WORLD, IERR)
      
      !iogni processore riceve il buffer e si costruisce il suo pezzo di matrice

      ind_loc = 1
      do jk =1 , jpk
          do jj =1 , jpj
                  do ji =1 , jpi
                          copy_in(ji,jj,jk)=buff_scatter(ind_loc)
                          ind_loc = ind_loc + 1
                  enddo
          enddo
      enddo

      do ji =1 , jpi
          do jj =1 , jpj
                  do jk =1 , jpk
                          matrix_scatter_local(jk,jj,ji)=real(copy_in(ji,jj,jk),8)
                  enddo
          enddo
      enddo

      END SUBROUTINE read_float

!-----------------------------------------------------------------------------------------------------------------------------------------

      SUBROUTINE read_float_2d(filename, var_to_store, shift,reading_proc,matrix_scatter_local_2d)
      USE calendar
      USE myalloc
      USE TIME_MANAGER
      USE ogstm_mpi_module

      USE MPI_GATHER_INFO
      USE MATRIX_VARS
      USE NODES_MODULE
      USE DTYPE_PROCS_STRING_MODULE

      character, INTENT(IN) :: filename*(*) ,var_to_store*(*)
      INTEGER, INTENT(IN) :: shift, reading_proc

      double precision,INTENT(OUT) :: matrix_scatter_local_2d(jpj,jpi)
      real, dimension(jpi,jpj) :: copy_in_2d

      real, dimension(jpiglo,jpjglo) :: buf_glo_2d

      INTEGER idrank, ierr, istart, jstart, status(MPI_STATUS_SIZE)
      INTEGER totistart, totiend
      INTEGER totjstart, totjend
      INTEGER jk,jj,ji,i,j,k
      INTEGER ind,ind_loc

      real, allocatable, dimension(:) :: buff_scatterv_tot_2d
      real, allocatable :: buff_scatter_2d(:)

      allocate(buff_scatter_2d (jpi* jpj))
      buff_scatter_2d = huge(buff_scatter_2d(1))
      ALLOCATE (buff_scatterv_tot_2d(cont_2d))
      buff_scatterv_tot_2d = huge(buff_scatterv_tot_2d(1))

      if(myrank.eq.0) then
              call readnc_float_global_2d(filename,var_to_store,buf_glo_2d,ingv_lon_shift)
      end if

      !!! il proc 0 deve trasformare la matrice 3d in un array monodimensionale, secondo gli indici dei processori

      if(myrank.eq.0) then
        ind=1
        DO idrank = 0,mpi_glcomm_size-1
                totistart = istart_a(idrank + 1)
                totiend   = istart_a(idrank + 1) + jpi_rec_a(idrank + 1) -1
                totjstart = jstart_a(idrank + 1)
                totjend   = jstart_a(idrank + 1) + jpj_rec_a(idrank + 1) -1

                do jj =totjstart,totjend
                        do ji =totistart,totiend
                                buff_scatterv_tot_2d(ind)=buf_glo_2d(ji,jj)
                                ind=ind+1
                        enddo
                enddo
        END DO
      end if


      !adesso il proc 0 manda i pezzi del buff scatterv tot a tutti i procs

      call mppsync()

      CALL MPI_SCATTERV(buff_scatterv_tot_2d, jprcv_count_2d,jpdispl_count_2d, MPI_REAL, buff_scatter_2d, sendcount_2d,MPI_REAL, 0, MPI_COMM_WORLD, IERR)

      !iogni processore riceve il buffer e si costruisce il suo pezzo di matrice

      ind_loc = 1
      do jj =1 , jpj
        do ji =1 , jpi
                copy_in_2d(ji,jj)=buff_scatter_2d(ind_loc)
                ind_loc = ind_loc + 1
        enddo
      enddo

      do ji =1 , jpi
          do jj =1 , jpj
                matrix_scatter_local_2d(jj,ji)=real(copy_in_2d(ji,jj),8)
          enddo
      enddo

      END SUBROUTINE read_float_2d

!----------------------------------------------------------------------------------------------------------

      SUBROUTINE read_double(filename, var_to_store,reading_proc,matrix_scatter_local)
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
      INTEGER, INTENT(IN) :: reading_proc

      double precision,INTENT(OUT) :: matrix_scatter_local(jpk,jpj,jpi)
      double precision, dimension(jpi,jpj,jpk) :: copy_in
      double precision, dimension(jpiglo,jpjglo,jpk) :: buf_glo

      INTEGER idrank, ierr, istart, jstart, status(MPI_STATUS_SIZE)
      INTEGER totistart, totiend
      INTEGER totjstart, totjend
      INTEGER jk,jj,ji,i,j,k
      INTEGER ind,ind_loc

      double precision, allocatable, dimension(:) :: buff_scatterv_tot
      double precision, allocatable :: buff_scatter(:)

      allocate(buff_scatter (jpi* jpj* jpk))
      buff_scatter = huge(buff_scatter(1))
      ALLOCATE (buff_scatterv_tot(total_dim_jpijk))
      buff_scatterv_tot = huge(buff_scatterv_tot(1))

      if(myrank.eq.0) then
              call readnc_double_global(filename,var_to_store,buf_glo)
      end if

      !!! il proc 0 deve trasformare la matrice 3d in un array monodimensionale, secondo gli indici dei processori

      if(myrank.eq.0) then
        ind=1
        DO idrank = 0,mpi_glcomm_size-1
                totistart = istart_a(idrank + 1)
                totiend   = istart_a(idrank + 1) + jpi_rec_a(idrank + 1) -1
                totjstart = jstart_a(idrank + 1)
                totjend   = jstart_a(idrank + 1) + jpj_rec_a(idrank + 1) -1

                do jk =1, jpk
                        do jj =totjstart,totjend
                                do ji =totistart,totiend
                                        buff_scatterv_tot(ind)=buf_glo(ji,jj,jk)
                                        ind=ind+1
                                enddo
                        enddo
                enddo
        END DO
      end if


      !adesso il proc 0 manda i pezzi del buff scatterv tot a tutti i procs

      call mppsync()

      CALL MPI_SCATTERV(buff_scatterv_tot, jprcv_count, jpdispl_count,MPI_DOUBLE, buff_scatter, sendcount,MPI_DOUBLE, 0, MPI_COMM_WORLD, IERR)

      !iogni processore riceve il buffer e si costruisce il suo pezzo di matrice

      ind_loc = 1
      do jk =1 , jpk
          do jj =1 , jpj
                  do ji =1 , jpi
                          copy_in(ji,jj,jk)=buff_scatter(ind_loc)
                          ind_loc = ind_loc + 1
                  enddo
          enddo
      enddo

      do ji =1 , jpi
          do jj =1 , jpj
                  do jk =1 , jpk
                          matrix_scatter_local(jk,jj,ji)=copy_in(ji,jj,jk)
                  enddo
          enddo
      enddo

      END SUBROUTINE read_double

!------------------------------------------------------------------------------------------------------------------

      SUBROUTINE read_double_2d(filename, var_to_store,reading_proc,matrix_scatter_local_2d)
      USE calendar
      USE myalloc
      USE TIME_MANAGER
      USE ogstm_mpi_module

      USE MPI_GATHER_INFO
      USE MATRIX_VARS
      USE NODES_MODULE
      USE DTYPE_PROCS_STRING_MODULE

      character, INTENT(IN) :: filename*(*) ,var_to_store*(*)
      INTEGER, INTENT(IN) :: reading_proc

      double precision,INTENT(OUT) :: matrix_scatter_local_2d(jpj,jpi)
      double precision, dimension(jpi,jpj) :: copy_in_2d

      double precision, dimension(jpiglo,jpjglo) :: buf_glo_2d

      INTEGER idrank, ierr, istart, jstart, status(MPI_STATUS_SIZE)
      INTEGER totistart, totiend
      INTEGER totjstart, totjend
      INTEGER jk,jj,ji,i,j,k
      INTEGER ind,ind_loc

      double precision, allocatable, dimension(:) :: buff_scatterv_tot_2d
      double precision, allocatable :: buff_scatter_2d(:)

      allocate(buff_scatter_2d (jpi* jpj))
      buff_scatter_2d = huge(buff_scatter_2d(1))
      ALLOCATE (buff_scatterv_tot_2d(cont_2d))
      buff_scatterv_tot_2d = huge(buff_scatterv_tot_2d(1))

      if(myrank.eq.0) then
              call readnc_double_global_2d(filename,var_to_store,buf_glo_2d)
      end if

      !!! il proc 0 deve trasformare la matrice 3d in un array monodimensionale, secondo gli indici dei processori

      if(myrank.eq.0) then
        ind=1
        DO idrank = 0,mpi_glcomm_size-1
                totistart = istart_a(idrank + 1)
                totiend   = istart_a(idrank + 1) + jpi_rec_a(idrank + 1) -1
                totjstart = jstart_a(idrank + 1)
                totjend   = jstart_a(idrank + 1) + jpj_rec_a(idrank + 1) -1

                do jj =totjstart,totjend
                        do ji =totistart,totiend
                                buff_scatterv_tot_2d(ind)=buf_glo_2d(ji,jj)
                                ind=ind+1
                        enddo
                enddo
        END DO
      end if


      !adesso il proc 0 manda i pezzi del buff scatterv tot a tutti i procs

      call mppsync()

      CALL MPI_SCATTERV(buff_scatterv_tot_2d,jprcv_count_2d,jpdispl_count_2d, MPI_DOUBLE, buff_scatter_2d, sendcount_2d,MPI_DOUBLE, 0, MPI_COMM_WORLD, IERR)

      !iogni processore riceve il buffer e si costruisce il suo pezzo di matrice

      ind_loc = 1
      do jj =1 , jpj
        do ji =1 , jpi
                copy_in_2d(ji,jj)=buff_scatter_2d(ind_loc)
                ind_loc = ind_loc + 1
        enddo
      enddo

      do ji =1 , jpi
          do jj =1 , jpj
                matrix_scatter_local_2d(jj,ji)=copy_in_2d(ji,jj)
          enddo
      enddo

      END SUBROUTINE read_double_2d

