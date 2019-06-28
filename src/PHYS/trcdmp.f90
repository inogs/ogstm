      SUBROUTINE trcdmp
!---------------------------------------------------------------------
!
!                       ROUTINE trcdmp
!                     ******************
!
!  Purpose :
!  --------
!    Compute (if asked) the passive tracer trend due to a newtonian
!    damping of the tracer field towards given data field and add it
!    to the general tracer trends.
!
!   Method :
!   --------
!    Default key          : empty routine, no damping trend
!    'key_trc_dmp' defined :
!       Newtonian damping towards tdta and sdta computed and add to
!       the general tracer trends:
!                     trc = ta + restotrc * (trcdta - trcb)
!       The trend is computed either throughout the water column
!       (nlmdmptrc=0) or in area of weak vertical mixing (nlmdmptrc=1) or
!       below the well mixed layer (nlmdmptrc=2)



       USE myalloc
       USE BC_mem
       use mpi
       IMPLICIT NONE

!----------------------------------------------------------------------
! local declarations
! ==================

      INTEGER :: jk,jj,ji, jn
      INTEGER :: tra_idx
      INTEGER :: jv
      INTEGER :: shift
      double precision :: ztra



      trcdmpparttime = MPI_WTIME() ! cronometer-start



! Boundary conditions for Atmosphere


       IF ( latmosph ) THEN
        DO jn=1, jn_atm
             tra_idx=tra_matrix_atm(jn)
          DO ji=1,jpi
        DO jj=1,jpj
        tra(1,jj,ji,tra_idx) = tra(1,jj,ji,tra_idx) + atm(jj,ji,jn)/e3t(1,jj,ji)
        ENDDO
        ENDDO
        ENDDO

       ENDIF




       trcdmpparttime = MPI_WTIME() - trcdmpparttime
       trcdmptottime  = trcdmptottime + trcdmpparttime



      END SUBROUTINE trcdmp
