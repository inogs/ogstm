PROGRAM OGSTM_MAIN
!
!
!                       PROGRAM MAIN
!                     ******************
!

      USE ogstm
      implicit none
      integer :: info,ierr
      write(*,*) "Starting ..."
      CALL mpi_init(ierr)

      !$OMP PARALLEL
      !$OMP MASTER
      CALL ogstm_launcher()
      !$OMP END MASTER
      !$OMP END PARALLEL

      CALL mpi_finalize(info)
      
END PROGRAM OGSTM_MAIN
