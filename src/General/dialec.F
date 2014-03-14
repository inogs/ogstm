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
      NAMELIST/frequency_avedump/ highfreq_tra_dia, freq_ave_phys


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
      clname = 'namelist.diagnostics'
      CALL ctlopn( numnam, clname, clold, clfor, clseq,ilseq, ilu, iused, numout, lwp, 1)


! 1. Read namlist files
! ---------------------

      DO ji=1,jptra_dia
        WRITE (dia(ji),1414) ji ! initalize diagnostic names for each slot
      END DO


      REWIND(numnam)
      READ(numnam,DIAGNOSTICS)
      REWIND(numnam)
      READ(numnam,frequency_avedump)

      jptra_dia_high= 0
      do ji =1, jptra_dia
          IF (highfreq_tra_dia(ji).eq.1) jptra_dia_high = jptra_dia_high + 1
      ENDDO
      allocate(highfreq_table_dia(jptra_dia_HIGH))

      jptra_dia_high = 0

      do ji =1, jptra_dia
          IF (highfreq_tra_dia(ji).eq.1) then
            jptra_dia_high = jptra_dia_high + 1
            highfreq_table_dia(jptra_dia_high) = ji
            WRITE(numout,*) dia(ji), ' belongs also to high freq group'
          ELSE
            WRITE(numout,*) dia(ji), ' belongs only to low freq group'
          ENDIF
      enddo

      if (freq_ave_phys.eq.2) then
           WRITE(numout,*) 'Forcings phys will be dumped following 2.aveTimes file'
           else
           if (freq_ave_phys.eq.2) then
             WRITE(numout,*) 'Forcings phys will be dumped following 1.aveTimes file'
           else
             WRITE(numout,*) 'Forcings phys will not be dumped'
           endif
      endif


      CLOSE( numnam )

      !CALL FLUSH(numout)


1414  FORMAT('D',I2.2)

      END SUBROUTINE dialec
