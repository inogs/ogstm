      PROGRAM OGSTM_MAIN
!
!                       PROGRAM MAIN
!                     ******************
!
!
      USE ogstm
      implicit none
      write(*,*) "Starting ..."

!$OMP PARALLEL
!$OMP MASTER
      CALL ogstm_launcher()
!$OMP END MASTER
!$OMP END PARALLEL

      END PROGRAM OGSTM_MAIN
