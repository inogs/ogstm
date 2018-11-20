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


! Boundary conditions fo Gib area

       !IF (Gsize .NE. 0) THEN
         !DO jn=1, jn_gib
          !tra_idx=tra_matrix_gib(jn)
          !DO jv=1, Gsize
             !ji = gib_ridxt(4,jv)
             !jj = gib_ridxt(3,jv)
             !jk = gib_ridxt(2,jv)
             !ztra = restotr(jk,jj,ji,tra_idx) * ( gib(jv,jn)-trb(jk,jj,ji,tra_idx) )
             !tra(jk,jj,ji,tra_idx) = tra(jk,jj,ji,tra_idx) + ztra
          !ENDDO
         !ENDDO
       !ENDIF


! Boundary conditions for rivers

       !if ( lrivers ) THEN
       !IF (Rsize .NE. 0) THEN
         !DO jn=1, jn_riv
          !tra_idx=tra_matrix_riv(jn)

          !DO jv=1, Rsize
             !ji = riv_ridxt(4,jv)
             !jj = riv_ridxt(3,jv)
            !tra(1,jj,ji,tra_idx) = tra(1,jj,ji,tra_idx) + riv(jv,jn)/e3t(1,jj,ji)
          !ENDDO
         !ENDDO
       !ENDIF
       !ENDIF


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
