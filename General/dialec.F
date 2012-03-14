! $Header: /cvsroot/opatm-bfm/opa_model/OPA/parlec.F,v 1.3 2010-03-05 09:56:44 cvsogs01 Exp $

      SUBROUTINE dialec
!---------------------------------------------------------------------
!
!                       ROUTINE parlec
!                     ******************
!
!  Purpose :
!  --------
!     Read and print options for the run (namelist)
!
!      the namelist file ( unit numnam ) :
!            &DIAGNOSTICS           : list of diagnostics
!
!   Output :
!   -------


       USE myalloc
       IMPLICIT NONE

! local declarations
! ==================


      CHARACTER*32 clname
      CHARACTER*21 clold, clfor, clseq, clnew, cldir, clunf
       integer ilseq,ji
      INTEGER ilu, iused(1,100)
      NAMELIST/DIAGNOSTICS/ dia



!----------------------------------------------------------------------
! statement functions
! ===================

#include "stafun.h"
! ... Open specifier
      clold = 'OLD'
      clnew = 'NEW'
      clfor = 'FORMATTED'
      clunf = 'UNFORMATTED'
      clseq = 'SEQUENTIAL'
      cldir = 'DIRECT'
      ilseq = 1 ! ... Sequential value

      DO ji = 1, 100
        iused(1,ji) = 0
      END DO


      IF(lwp) THEN
          WRITE(numout,*) ' '
          WRITE(numout,*) ' '
          WRITE(numout,*) ' routine dialec'
          WRITE(numout,*) ' **************'
          WRITE(numout,*) ' '
          WRITE(numout,*) ' initialization of run'
          WRITE(numout,*) ' '
          WRITE(numout,*) ' '
          WRITE(numout,*) ' namelist.diagnostics'
          WRITE(numout,*) ' ********'
          WRITE(numout,*) ' '
      ENDIF

      numnam = 2004
!      REWIND( numnam )
      clname = 'namelist.diagnostics'
      CALL ctlopn( numnam, clname, clold, clfor, clseq,ilseq, ilu, iused, numout, lwp, 1)


! 1. Read namlist files
! ---------------------

      DO ji=1,jptra_dia
        WRITE (dia(ji),'("TR",I1)') ji
      END DO


         REWIND(numnat)
         READ(numnat,dia)



      CLOSE( numnam )
      CALL FLUSH_(numout)
      END SUBROUTINE dialec
