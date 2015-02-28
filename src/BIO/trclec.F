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
!            &natrtd           : dynamical trends
!            &natadd           : additional arrays for diagnostics
!
!   OUTPUT :
!   ------



       USE myalloc
       IMPLICIT NONE


!----------------------------------------------------------------------
! local declarations
! ==================

      INTEGER ji
      INTEGER ilu,iused(1,100),ilseq
      CHARACTER(LEN=21) clold,clfor,clseq,clnew,cldir,clunf
      CHARACTER(LEN=32) clname
!----------------------------------------------------------------------
! statement functions
! ===================



C 0. initializations

      namelist /NATTRC/       ctrcnm, ctrcun, ctrmax,ctr_hf ! BUONO
      namelist /NATTRC_DIAG/   dianm,  diaun,         diahf ! BUONO


#    if defined key_trc_diatrd
      namelist /natrtd/ luttrd,nwritetrd                      ! dynamical trends
#    endif

#    if defined key_trc_diaadd
      namelist /natadd/ ctrc3d,ctrc3l,ctrc2d,ctrc2l, ctrc3u, ctrc2u,nwriteadd    !additional diagnostics
#    endif



! OPEN specifier

      clold = 'old'
      clnew = 'new'
      clfor = 'formatted'
      clunf = 'unformatted'
      clseq = 'SEQUENTIAL'
      cldir = 'direct'
      ilseq =  1! SEQUENTIAL value
      ilu =    0! initialize the number of LOGICAL UNIT used

      IF(lwp) THEN
          WRITE(numout,*) ' '
          WRITE(numout,*) ' ROUTINE trclec'
          WRITE(numout,*) ' **************'
          WRITE(numout,*) ' '
          WRITE(numout,*) ' namelist for passive tracers'
          WRITE(numout,*) ' ****************************'
          WRITE(numout,*) ' '
      ENDIF

      numnat=80

      clname='namelist.passivetrc'
      CALL ctlopn(numnat,clname,clold,clfor,clseq,ilseq,ilu,iused,numout,lwp,1)



!      ***********  namelist nattrc STATE VARIABLES:



      REWIND(numnat)
      READ(numnat,nattrc)

      IF(lwp) THEN
          WRITE(numout,*) ' '
          WRITE(numout,*) 'nattrc'
          WRITE(numout,*) ' '
!          WRITE(numout,*)' frequency of outputs for passive tracers nwritetrc = ',nwritetrc
!          WRITE(numout,*) ' restart LOGICAL for passive tr. lrsttr = ',lrsttr
!          WRITE(numout,*) ' control of time step for p. tr. nrsttr = ',nrsttr
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






#    if defined key_trc_diatrd

! 1.2 namelist nattrd : passive tracers dynamical trends

      nwritetrd = 10

! default : no dynamical trend recording

      DO ji=1,jptra
        luttrd(ji) = .FALSE.
      END DO 

      REWIND(numnat)
      READ(numnat,natrtd)

      IF(lwp) THEN
          WRITE(numout,*) 'natrtd'
          WRITE(numout,*) ' '
          WRITE(numout,*)' frequency of outputs for dynamical trends nwritetrd = ',nwritetrd
          DO ji=1,jptra
            WRITE(numout,*)' keep dynamical trends for tracer number :',ji,luttrd(ji)
          END DO 
      ENDIF


! 1.3 namelist natadd : passive tracers diagnostics



      nwriteadd = 10

! default value for 3D output arrays : short and long name, units

      DO ji=1,jpdia3d
        WRITE (ctrc3d(ji),'("3D_",I1)') ji
        WRITE (ctrc3l(ji),'("3D DIAGNOSTIC NUMBER ",I2)') ji
        ctrc3u(ji)=' '
      END DO 


! default value for 2D output arrays : short and long name, units

      DO ji=1,jpdia2d
        WRITE (ctrc2d(ji),'("2D_",I1)') ji
        WRITE (ctrc2l(ji),'("2D DIAGNOSTIC NUMBER ",I2)') ji
        ctrc2u(ji)=' '
      END DO 

      REWIND(numnat)
      READ(numnat,natadd)

      IF(lwp) THEN
          WRITE(numout,*) ' natadd'
          WRITE(numout,*) ' '
          WRITE(numout,*)' frequency of outputs for additional arrays nwriteadd = ',nwriteadd
          DO ji=1,jpdia3d
            WRITE(numout,*)'name of 3d output field number :',ji,' : ',ctrc3d(ji)
            WRITE(numout,*) ctrc3l(ji)  
            WRITE(numout,*) ' in unit = ',ctrc3u(ji)
          END DO 
          WRITE(numout,*) ' '
          DO ji=1,jpdia2d
            WRITE(numout,*) 'name of 2d output field number :',ji,' : ',ctrc2d(ji)
            WRITE(numout,*) ctrc2l(ji)  
            WRITE(numout,*) ' in unit = ',ctrc2u(ji)
          END DO 
          WRITE(numout,*) ' '
      ENDIF
#    endif


      CLOSE(numnat)
      END SUBROUTINE trclec
