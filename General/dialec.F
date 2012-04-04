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


      CHARACTER(LEN=32) clname
      CHARACTER(LEN=21) clold, clfor, clseq, clnew, cldir, clunf
       integer ilseq,ji
      INTEGER ilu, iused(1,100)
      NAMELIST/DIAGNOSTICS/ dia



!----------------------------------------------------------------------
! statement functions
! ===================

! ... Open specifier
      clold = 'OLD'
      clnew = 'NEW'
      clfor = 'FORMATTED'
      clunf = 'UNFORMATTED'
      clseq = 'SEQUENTIAL'
      cldir = 'DIRECT'
      ilseq = 1 ! ... Sequential value
      ilu = 0
      ilu=0

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
      REWIND( numnam )
      clname = 'namelist.diagnostics'
      CALL ctlopn( numnam, clname, clold, clfor, clseq,ilseq, ilu, iused, numout, lwp, 1)


! 1. Read namlist files
! ---------------------

      DO ji=1,jptra_dia
        WRITE (dia(ji),1414) ji ! initalize diagnostic names for each slot
      END DO


      REWIND(numnam)
      READ(numnam,DIAGNOSTICS)



      CLOSE( numnam )

      CALL FLUSH(numout)


1414  FORMAT('D',I2.2)

      END SUBROUTINE dialec
