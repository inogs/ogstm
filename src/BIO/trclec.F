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
!            &natnum           : numerical schemes
!            &natrtd           : dynamical trends
!            &natadd           : additional arrays for diagnostics
!
!   OUTPUT :
!   ------
!      COMMON                  : (#ifdef "key_passivetrc")
!            /cittrc/          : tracer names
!            /cot3ad/          : options for smolarkiewicz scheme
!            /cit3ad/          : other options for smolarkiewicz scheme
!            /cottrd/          : passive tracer trends (#ifef
!                                 "key_trc_diatrd")
!
!   WORKSPACE :                : no


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

      namelist /nattrc/ nwritetrc,lrsttr,nrsttr, ctrcnm,ctrcnl,ctrcun,lutini
      namelist /natnum/ rsc,rtrn,ncor,ndttrc,crosster,lhdf,ahtrb0,trcrat,ahtrc0,aeivtr0,vsed,photop,atlantic_bfm,bottom_flux

#    if defined key_trc_diatrd
      namelist /natrtd/ luttrd,nwritetrd                      ! dynamical trends
#    endif

#    if defined key_trc_diaadd
      namelist /natadd/ ctrc3d,ctrc3l,ctrc2d,ctrc2l, ctrc3u, ctrc2u,nwriteadd    !additional diagnostics
#    endif
#    if defined key_trc_dmp
!      namelist/natdmp/nmldmptr
#    endif
      namelist /MAXVALUES/isCheckLOG, maxvalue
      namelist /frequency_avedump/ highfreq_tra

! OPEN specifier

      clold = 'old'
      clnew = 'new'
      clfor = 'formatted'
      clunf = 'unformatted'
      clseq = 'SEQUENTIAL'
      cldir = 'direct'

! SEQUENTIAL value

      ilseq=1

! initialize the number of LOGICAL UNIT used

      ilu=0

      IF(lwp) THEN
          WRITE(numout,*) ' '
          WRITE(numout,*) ' ROUTINE trclec'
          WRITE(numout,*) ' **************'
          WRITE(numout,*) ' '
          WRITE(numout,*) ' namelist for passive tracer'
          WRITE(numout,*) ' ***************************'
          WRITE(numout,*) ' '
      ENDIF

      numnat=80

      clname='namelist.passivetrc'
      CALL ctlopn(numnat,clname,clold,clfor,clseq,ilseq,ilu,iused,numout,lwp,1)


! 1., 2. & 3. initialization with namelist files
! ----------------------------------------------

! 1.0 namelist nattrc :

      nwritetrc = 10
      lrsttr    =.FALSE.
      nrsttr    = 0

      DO ji=1,jptra

!        WRITE (ctrcnm(ji),'("TR_",I1)') ji
        WRITE (ctrcnl(ji),'("TRACER NUMBER ",I1)') ji
        ctrcun(ji)='mmole/m3'
        lutini(ji)=.FALSE. 
      END DO


         REWIND(numnat)
         READ(numnat,nattrc)

      IF(lwp) THEN
          WRITE(numout,*) ' '
          WRITE(numout,*) 'nattrc'
          WRITE(numout,*) ' '
          WRITE(numout,*)' frequency of outputs for passive tracers nwritetrc = ',nwritetrc
          WRITE(numout,*) ' restart LOGICAL for passive tr. lrsttr = ',lrsttr
          WRITE(numout,*) ' control of time step for p. tr. nrsttr = ',nrsttr
          DO ji=1,jptra
            WRITE(numout,*) ' tracer nb: ',ji,' name = ',ctrcnm(ji),ctrcnl(ji)
            WRITE(numout,*) ' in unit = ',ctrcun(ji)
            WRITE(numout,*) ' initial value in FILE : ',lutini(ji) 
            WRITE(numout,*) ' '
          END DO
          WRITE(numout,*) ' '
      ENDIF

! 1.1 namelist natnum :

      rsc         = 1.
      rtrn        = 1.e-15
      ncor        = 1
      ndttrc      = 4
      crosster    = .FALSE.
      lhdf        = .TRUE.
      ahtrb0      = 0.
      trcrat      = 1.
      ahtrc0      = aht0
      aeivtr0     = aeiv0
      vsed        = 3.0
      photop      = .FALSE.
      atlantic_bfm= .FALSE.
      bottom_flux = 0.

      REWIND(numnat)
      READ(numnat,natnum)

      IF(lwp) THEN
          WRITE(numout,*) ' '
          WRITE(numout,*) 'natnum'
          WRITE(numout,*) ' '
          WRITE(numout,*) ' rsc tuning coefficient                            = ' , rsc
          WRITE(numout,*) ' rtrn truncation value                             = ' , rtrn
          WRITE(numout,*) ' ncor number of corrective phase                   = ', ncor
          WRITE(numout,*) ' ndttrc time step freq. for pass. trac.            = ', ndttrc
          WRITE(numout,*) ' crosster  computes or not crossterms              = ', crosster
          WRITE(numout,*) ' lhdf  calls or not trchdf                         = ', lhdf
          WRITE(numout,*) ' background diffusivity for passive tr             = ', ahtrb0
          WRITE(numout,*) ' ratio betweeen passive and active tr diffusion coeff= ', trcrat
          WRITE(numout,*) ' horizontal eddy diffus. for passive tr            = ', ahtrc0
          WRITE(numout,*) ' eddy induced veloc. coef. for passive tr          = ', aeivtr0
          WRITE(numout,*) ' detritus sedimentation speed   vsed               =', vsed/86400
          WRITE(numout,*) ' photoperiod scaling photop                        =', photop
          WRITE(numout,*) ' activation of bfm in atlantic buffer              =', atlantic_bfm
          WRITE(numout,*) ' bottom flux [0,1], 0 -> no flux, 1 -> total flux  =', bottom_flux
      ENDIF

      IF (vsed .LT. 0.) THEN
          write (*,*) 'vsed must be greated than 0 instead it is:', vsed/86400
          STOP
      ENDIF

      IF ((bottom_flux .LT. 0.) .OR. (bottom_flux .GT. 1.)) THEN
          write (*,*) 'bottom flux must be in [0,1] instead it is:', bottom_flux
          STOP
      ENDIF
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


      REWIND(numnat)
      READ(numnat,maxvalues)
      REWIND(numnat)
      READ(numnat,frequency_avedump)
      CLOSE(numnat)

      jptra_high = 0
      do ji =1,jptra
          if (highfreq_tra(ji).eq.1) jptra_high = jptra_high + 1
      enddo
      allocate(highfreq_table(jptra_HIGH))

      jptra_high = 0
      do ji =1,jptra
          if (highfreq_tra(ji).eq.1) then
              jptra_high = jptra_high + 1
              highfreq_table(jptra_high) = ji
              WRITE(numout,*) ctrcnm(ji), ' belongs also to high freq group'
           else
               WRITE(numout,*) ctrcnm(ji), ' belongs only to low freq group'
           endif
      enddo

      END SUBROUTINE trclec
