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
!   INPUT :
!   -----
!      the namelist FILE ( UNIT numnat ) :
!            &nattrc           : general
!            &nattrc_diag      : general

       USE myalloc
       IMPLICIT NONE

!----------------------------------------------------------------------
! local declarations
! ==================

      INTEGER ji

!----------------------------------------------------------------------
! statement functions
! ===================

      namelist /NATTRC/       ctrcnm, ctrcun, ctrmax,ctr_hf
      namelist /NATTRC_DIAG/   dianm,  diaun,         diahf
      namelist /NATTRC_DIAG_2d/   dianm_2d,  diaun_2d,  diahf_2d

      IF(lwp) THEN
          WRITE(numout,*) ' '
          WRITE(numout,*) ' ROUTINE trclec'
          WRITE(numout,*) ' **************'
          WRITE(numout,*) ' '
          WRITE(numout,*) ' namelist for passive tracers'
          WRITE(numout,*) ' ****************************'
          WRITE(numout,*) ' '
      ENDIF

      OPEN(unit=numnat, file='namelist.passivetrc', status= 'OLD') ! 'FORMATTED', 'SEQUENTIAL')




!      ***********  namelist nattrc STATE VARIABLES:



      REWIND(numnat)
      READ(numnat,nattrc)

      IF(lwp) THEN
          WRITE(numout,*) ' '
          WRITE(numout,*) 'nattrc'
          WRITE(numout,*) ' '
          DO ji=1,jptra
            WRITE(numout,250) 'tracer nb: ',ji,' name = ',ctrcnm(ji), ' in unit = ',ctrcun(ji), ' max ',
     +       ctrmax(ji), ' highFreq ', ctr_hf(ji)
250   FORMAT (A,I3, A8,A15,A,A,A,ES20.7,A,I1) 
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
              if (lwp) WRITE(numout,*) ctrcnm(ji), ' belongs also to high freq group'
           else
               if (lwp) WRITE(numout,*) ctrcnm(ji), ' belongs only to low freq group'
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




!      ***********  namelist nattrc_diag_2d DIAGNOSTIC VARIABLES surface or bottom:

      REWIND(numnat)
      READ(numnat,nattrc_diag_2d)

      jptra_dia2d_high= 0
      do ji =1, jptra_dia_2d
          IF (diahf_2d(ji).eq.1) jptra_dia2d_high = jptra_dia2d_high + 1
      ENDDO

      if (lwp) write(*,*) 'High freq diagnostics number 2d:', jptra_dia2d_HIGH
      allocate(highfreq_table_dia2d(jptra_dia2d_HIGH))

      jptra_dia2d_high = 0

      do ji =1, jptra_dia_2d
          IF (diahf_2d(ji).eq.1) then
            jptra_dia2d_high = jptra_dia2d_high + 1
            highfreq_table_dia2d(jptra_dia2d_high) = ji
            if (lwp) WRITE(numout,*) dianm_2d(ji), ' belongs also to high freq group'
          ELSE
            if (lwp) WRITE(numout,*) dianm_2d(ji), ' belongs only to low freq group'
          ENDIF
      enddo




      CLOSE(numnat)


      IF(lwp) THEN
          WRITE(numout,*) ' '
          WRITE(numout,*) ' *** number of passive tracer jptra = ',jptra
          WRITE(numout,*) ' '
      ENDIF

      END SUBROUTINE trclec
