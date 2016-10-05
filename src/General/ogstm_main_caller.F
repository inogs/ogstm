      PROGRAM OGSTM_MAIN
!
!                       PROGRAM MAIN
!                     ******************
!
!
!
      write(*,*) "Starting ..."

!$OMP PARALLEL
!$OMP MASTER
      CALL ogstm()
!$OMP END MASTER
!$OMP END PARALLEL

      END PROGRAM OGSTM_MAIN
