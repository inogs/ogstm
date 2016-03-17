      SUBROUTINE benlec
!---------------------------------------------------------------------
!
!                       ROUTINE benlec
!                     ******************
!
!  PURPOSE :
!  ---------
!     READ and PRINT options for the passive tracer run (namelist)
!
!   INPUT :
!   -----
!      the namelist FILE ( UNIT numben ) :
!            &nattrc           : general
!            &nattrc_diag      : general

       USE myalloc
       USE BIOBEN_mem
       IMPLICIT NONE

!----------------------------------------------------------------------
! local declarations
! ==================

      INTEGER ji

!----------------------------------------------------------------------
! statement functions
! ===================

      namelist /NATBEN/       cbennm, cbenun, cbemax,cbe_hf
      namelist /NATBEN_DIAG/   dibnm,  dibun,         dibhf
      namelist /NATBEN_DIAG_2d/   dibnm_2d,  dibun_2d,  dibhf_2d

      IF(lwp) THEN
          WRITE(numout,*) ' '
          WRITE(numout,*) ' ROUTINE benlec'
          WRITE(numout,*) ' **************'
          WRITE(numout,*) ' '
          WRITE(numout,*) ' namelist for benthic model'
          WRITE(numout,*) ' ****************************'
          WRITE(numout,*) ' '
      ENDIF

      OPEN(unit=numben, file='namelist.ben', status= 'OLD') ! 'FORMATTED', 'SEQUENTIAL')




!      ***********  namelist natben STATE VARIABLES:



      REWIND(numben)
      READ(numben,natben)

      IF(lwp) THEN
          WRITE(numout,*) ' '
          WRITE(numout,*) 'natben'
          WRITE(numout,*) ' '
          DO ji=1,jptra_b
            WRITE(numout,250) 'tracer nb: ',ji,' name = ',cbennm(ji), ' in unit = ',cbenun(ji), ' max ',
     +       cbemax(ji), ' highFreq ', cbe_hf(ji)
250   FORMAT (A,I3, A8,A15,A,A,A,ES20.7,A,I1) 
            WRITE(numout,*) ' '
          END DO
          WRITE(numout,*) ' '
      ENDIF

      jptra_high_ben = 0
      do ji =1,jptra_b
          if (cbe_hf(ji).eq.1) jptra_high_ben = jptra_high_ben + 1
      enddo
      allocate(highfreq_table_ben(jptra_high_ben)) ; highfreq_table_ben = huge(highfreq_table_ben(1))

      jptra_high_ben = 0
      do ji =1,jptra_b
          if (cbe_hf(ji).eq.1) then
              jptra_high_ben = jptra_high_ben + 1
              highfreq_table_ben(jptra_high_ben) = ji
              if (lwp) WRITE(numout,*) cbennm(ji), ' belongs also to high freq group'
           else
               if (lwp) WRITE(numout,*) cbennm(ji), ' belongs only to low freq group'
           endif
      enddo


!      ***********  namelist natben_diag DIAGNOSTIC VARIABLES:

      REWIND(numben)
      READ(numben,natben_diag)

      jptra_dia_high_ben= 0
      do ji =1, jptra_dia_b
          IF (dibhf(ji).eq.1) jptra_dia_high_ben = jptra_dia_high_ben + 1
      ENDDO

      if (lwp) write(*,*) 'High freq diagnostics number :',  jptra_dia_high_ben
      allocate(highfreq_table_dia_ben(jptra_dia_high_ben))

      jptra_dia_high_ben = 0

      do ji =1, jptra_dia_b
          IF (dibhf(ji).eq.1) then
            jptra_dia_high_ben = jptra_dia_high_ben + 1
            highfreq_table_dia_ben(jptra_dia_high_ben) = ji
            if (lwp) WRITE(numout,*) dibnm(ji), ' belongs also to high freq group'
          ELSE
            if (lwp) WRITE(numout,*) dibnm(ji), ' belongs only to low freq group'
          ENDIF
      enddo




!      ***********  namelist natben_diag_2d DIAGNOSTIC VARIABLES surface or bottom:

      REWIND(numben)
      READ(numben,natben_diag_2d)

      jptra_dia2d_high_ben= 0
      do ji =1, jptra_dia_b_2d
          IF (dibhf_2d(ji).eq.1) jptra_dia2d_high_ben = jptra_dia2d_high_ben + 1
      ENDDO

      if (lwp) write(*,*) 'High freq diagnostics number 2d:', jptra_dia2d_high_ben
      allocate(highfreq_table_dia2d_ben(jptra_dia2d_high_ben))

      jptra_dia2d_high_ben = 0

      do ji =1, jptra_dia_b_2d
          IF (dibhf_2d(ji).eq.1) then
            jptra_dia2d_high_ben = jptra_dia2d_high_ben + 1
            highfreq_table_dia2d_ben(jptra_dia2d_high_ben) = ji
            if (lwp) WRITE(numout,*) dibnm_2d(ji), ' belongs also to high freq group'
          ELSE
            if (lwp) WRITE(numout,*) dibnm_2d(ji), ' belongs only to low freq group'
          ENDIF
      enddo




      CLOSE(numben)


      IF(lwp) THEN
          WRITE(numout,*) ' '
          WRITE(numout,*) ' *** number of benthic variables jptra_b = ',jptra_b
          WRITE(numout,*) ' '
      ENDIF

      END SUBROUTINE benlec
