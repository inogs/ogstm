C $Id: trcdmp.F,v 1.6 2010-05-26 10:38:12 cvsogs01 Exp $

      SUBROUTINE trcdmp
!---------------------------------------------------------------------
!
!                       ROUTINE trcdmp
!                     ******************
!
!  Purpose :
!  --------
!	Compute (if asked) the passive tracer trend due to a newtonian
!	damping of the tracer field towards given data field and add it
!	to the general tracer trends.
!
!   Method :
!   --------
!	Default key          : empty routine, no damping trend
!	'key_trc_dmp' defined :
!	   Newtonian damping towards tdta and sdta computed and add to
!	   the general tracer trends:
!                     trc = ta + restotrc * (trcdta - trcb)
!	   The trend is computed either throughout the water column
!	   (nlmdmptrc=0) or in area of weak vertical mixing (nlmdmptrc=1) or
!	   below the well mixed layer (nlmdmptrc=2)
!
!      'key_diatrdtrc' defined: the tracers damping trends are saved.
!
!      Macro-tasked on tracer (jn-loop)
!
!   Input :
!   ------

!      common
!            /cottrc/          : passive tracers fields (before,now,after)
!            /cotrda/          : passive tracers DATA READ and at given time_step
!            /comtsk/          : multitasking
!
!   Output :
!   -------
!      common
!	      /comaft/ ta,sa	: general tracer trend increased by
!				  the tracer newtonian damping trends
!	      /comtra/ ttrd,strd  tracers damping trends
!                                (IF 'key-diatrdtra' is defined)
!



       USE myalloc
       USE myalloc_mpp
       USE BC_mem
       IMPLICIT NONE

!----------------------------------------------------------------------
! local declarations
! ==================

      INTEGER ji, jj, jk, jn
!      INTEGER jn_nudg
      INTEGER tra_idx
      INTEGER nmldmp_gib,nmldmp_riv,nmldmp_atm
      INTEGER(4) jv
      REAL(8)  ztra
!      REAL(8) cvmgm, cvmgp, ztest
!----------------------------------------------------------------------
! statement functions
! ===================

#include "stafun.h"

! tracer loop parallelized (macrotasking)
! =======================================


       trcdmpparttime = MPI_WTIME() ! F79 cronometer-start


!      jn_nudg = 4
!      tra_matrix(1) = ppO2o
!      tra_matrix(2) = ppN1p
!      tra_matrix(3) = ppN3n
!      tra_matrix(4) = ppN5s

      nmldmp_gib = 0
      nmldmp_riv = 1
      nmldmp_atm = 1


! Boundary conditions fo Gib area
       IF (Gsize .NE. 0) THEN
         DO jn=1, jn_gib
          tra_idx=tra_matrix_gib(jn)
          DO jv=1, Gsize
             ji = gib_ridxt(jv,2) 
             jj = gib_ridxt(jv,3)
             jk = gib_ridxt(jv,4) 
           ztra = 0
           IF ( nmldmp_gib.EQ.0 ) THEN
             ztra = restotr(ji,jj,jk,tra_idx) * ( gib(jv,jn)-trb(ji,jj,jk,tra_idx) )
           ENDIF
           IF ( nmldmp_gib.EQ.1 ) THEN
             ztra = gib(jv,jn)
           ENDIF
             tra(ji,jj,jk,tra_idx) = tra(ji,jj,jk,tra_idx) + ztra
          ENDDO 
         ENDDO 
!        write(103+rank,*) 'NEW - > restotr', restotr(20,62,1,tra_matrix(1))
!        write(103+rank,*) 'NEW - > gib    ', gib(1,1)
!        write(103+rank,*) 'NEW - > trb    ', trb(20,62,1,tra_matrix(1))
!        write(103+rank,*) 'NEW - > ztra   ' , restotr(20,62,1,tra_matrix(1))* ( gib(1,1)-trb(20,62,1,tra_matrix(1)) )
!        CALL flush_(103+rank)
 
       ENDIF
! Boundary conditions for rivers
       IF (Rsize .NE. 0) THEN
         DO jn=1, jn_riv
          tra_idx=tra_matrix_riv(jn)
          DO jv=1, Rsize
             ji = riv_ridxt(jv,2)
             jj = riv_ridxt(jv,3)
             jk = riv_ridxt(jv,4)
           ztra = 0
           IF ( nmldmp_riv.EQ.0 ) THEN
             ztra = restotr(ji,jj,jk,tra_idx)* ( riv(jv,jn)-trb(ji,jj,jk,tra_idx) )
           ENDIF
           IF ( nmldmp_riv.EQ.1 ) THEN
             ztra = riv(jv,jn)
           ENDIF
             tra(ji,jj,jk,tra_idx) = tra(ji,jj,jk,tra_idx) + ztra
          ENDDO
         ENDDO
!        write(103+rank,*) 'OLD - >', restotr(20,62,1,tra_matrix(1)) * ( trdta(20,62,1,tra_matrix(1))-trb(20,62,1,tra_matrix(1)) )
!        write(103+rank,*) 'NEW - >', restotr(20,62,1,tra_matrix(1))* ( riv(1,1)-trb(20,62,1,tra_matrix(1)) )
!        CALL flush_(103+rank)

       ENDIF

! Boundary conditions for Atmosphere
       IF (Asize .NE. 0) THEN
         DO jn=1, jn_atm
          tra_idx=tra_matrix_atm(jn)
          DO jv=1, Asize
             ji = atm_ridxt(jv,2)
             jj = atm_ridxt(jv,3)
             jk = atm_ridxt(jv,4)
           ztra = 0
           IF ( nmldmp_atm.EQ.0 ) THEN
             ztra = restotr(ji,jj,jk,tra_idx)* ( atm(jv,jn)-trb(ji,jj,jk,tra_idx) )
           ENDIF
           IF ( nmldmp_atm.EQ.1 ) THEN
             ztra = atm(jv,jn)
           ENDIF
             tra(ji,jj,jk,tra_idx) = tra(ji,jj,jk,tra_idx) + ztra
          ENDDO
         ENDDO

       ENDIF




       trcdmpparttime = MPI_WTIME() - trcdmpparttime ! F79 cronometer-stop
       trcdmptottime = trcdmptottime + trcdmpparttime



      END SUBROUTINE trcdmp
