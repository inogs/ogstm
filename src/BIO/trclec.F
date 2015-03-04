      SUBROUTINE trclec
!---------------------------------------------------------------------
!
!                       ROUTINE trclec
!                     ******************
!
!  PURPOSE :
!  ---------
!     READ and PRINT options for the passive tracer run (namelist)
!
!   METHOD :                   : no
!   -------
!
!   INPUT :
!   -----
!      the namelist FILE ( UNIT numnat ) :
!            &nattrc           : general
!            &nattrc_diag      : general
!
!   OUTPUT :
!   ------



       USE myalloc
       IMPLICIT NONE


!----------------------------------------------------------------------
! local declarations
! ==================

      INTEGER ji
      INTEGER iused(1,100), ilseq,ilu
      CHARACTER(LEN=21) clold,clfor,clseq,clnew,cldir,clunf
      CHARACTER(LEN=32) clname
!----------------------------------------------------------------------
! statement functions
! ===================

      namelist /NATTRC/       ctrcnm, ctrcun, ctrmax,ctr_hf
      namelist /NATTRC_DIAG/   dianm,  diaun,         diahf



! OPEN specifier
      ilseq = 1
      ilu   = 0
      clold = 'old'
      clnew = 'new'
      clfor = 'formatted'
      clunf = 'unformatted'
      clseq = 'SEQUENTIAL'
      cldir = 'direct'

      IF(lwp) THEN
          WRITE(numout,*) ' '
          WRITE(numout,*) ' ROUTINE trclec'
          WRITE(numout,*) ' **************'
          WRITE(numout,*) ' '
          WRITE(numout,*) ' namelist for passive tracers'
          WRITE(numout,*) ' ****************************'
          WRITE(numout,*) ' '
      ENDIF

      clname='namelist.passivetrc'
      CALL ctlopn(numnat,clname,clold,clfor,clseq,ilseq,ilu,iused,numout,lwp,1)



!      ***********  namelist nattrc STATE VARIABLES:



      REWIND(numnat)
      READ(numnat,nattrc)

      IF(lwp) THEN
          WRITE(numout,*) ' '
          WRITE(numout,*) 'nattrc'
          WRITE(numout,*) ' '
          DO ji=1,jptra
            WRITE(numout,*) ' tracer nb: ',ji,' name = ',ctrcnm(ji), ' in unit = ',ctrcun(ji), ' max ',
     +       ctrmax(ji), 'highFreq ', ctr_hf(ji)
            WRITE(numout,*) ' '
          END DO
          WRITE(numout,*) ' '
      ENDIF

      jptra_high = 0
      do ji =1,jptra
          if (ctr_hf(ji).eq.1) jptra_high = jptra_high + 1
      enddo
      allocate(highfreq_table(jptra_HIGH)) ; highfreq_table = huge(highfreq_table(1))

      jptra_high = 0
      do ji =1,jptra
          if (ctr_hf(ji).eq.1) then
              jptra_high = jptra_high + 1
              highfreq_table(jptra_high) = ji
              WRITE(numout,*) ctrcnm(ji), ' belongs also to high freq group'
           else
               WRITE(numout,*) ctrcnm(ji), ' belongs only to low freq group'
           endif
      enddo


!      ***********  namelist nattrc_diag DIAGNOSTIC VARIABLES:

      REWIND(numnat)
      READ(numnat,nattrc_diag)

      jptra_dia_high= 0
      do ji =1, jptra_dia
          IF (diahf(ji).eq.1) jptra_dia_high = jptra_dia_high + 1
      ENDDO

      if (lwp) write(*,*) 'High freq diagnostics number :', jptra_dia_HIGH
      allocate(highfreq_table_dia(jptra_dia_HIGH))

      jptra_dia_high = 0

      do ji =1, jptra_dia
          IF (diahf(ji).eq.1) then
            jptra_dia_high = jptra_dia_high + 1
            highfreq_table_dia(jptra_dia_high) = ji
            if (lwp) WRITE(numout,*) dianm(ji), ' belongs also to high freq group'
          ELSE
            if (lwp) WRITE(numout,*) dianm(ji), ' belongs only to low freq group'
          ENDIF
      enddo



      CLOSE(numnat)
      END SUBROUTINE trclec
