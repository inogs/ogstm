PROGRAM OGSTM_MAIN
!
!
!                       PROGRAM MAIN
!                     ******************
!

      USE ogstm
#ifdef _OPENACC
      use openacc
#endif
      implicit none
      integer :: info, ierr

#ifdef _OPENACC
      integer :: num_gpus, gpu_rank, mpi_rank, acc_device_type
#endif
      
      write(*,*) "Starting ..."
      CALL mpi_init(ierr)

#ifdef _OPENACC
      acc_device_type = acc_get_device_type()
      num_gpus = acc_get_num_devices(acc_device_type)
      call MPI_COMM_RANK(MPI_COMM_WORLD, mpi_rank, ierr)
      call acc_set_device_num(mpi_rank, acc_device_type)
      gpu_rank =  acc_get_device_num(acc_device_type)
      print *, "My MPI rank: ", mpi_rank, "My GPU rank: ", gpu_rank, "out of ", num_gpus
#endif

      !$OMP PARALLEL
      !$OMP MASTER
      CALL ogstm_launcher()
      !$OMP END MASTER
      !$OMP END PARALLEL

      CALL mpi_finalize(info)
      
END PROGRAM OGSTM_MAIN
