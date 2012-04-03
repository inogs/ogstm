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
      CHARACTER*21 clold,clfor,clseq,clnew,cldir,clunf
      CHARACTER*32 clname
!----------------------------------------------------------------------
! statement functions
! ===================



C 0. initializations

      namelist /nattrc/ nwritetrc,lrsttr,nrsttr, ctrcnm,ctrcnl,ctrcun,lutini
      namelist /natnum/ rsc,rtrn,ncor,ndttrc,crosster,lhdf,ahtrb0,trcrat,ahtrc0,aeivtr0

#    if defined key_trc_diatrd
      namelist /natrtd/ luttrd,nwritetrd                      ! dynamical trends
#    endif

#    if defined key_trc_diaadd
      namelist /natadd/ ctrc3d,ctrc3l,ctrc2d,ctrc2l, ctrc3u, ctrc2u,nwriteadd    !additional diagnostics
#    endif
#    if defined key_trc_dmp
!      namelist/natdmp/nmldmptr
#    endif


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
      REWIND (numnat)

      clname='namelist.passivetrc'
      CALL ctlopn(numnat,clname,clold,clfor,clseq,ilseq,ilu,iused,numout,lwp,1)


! 1., 2. & 3. initialization with namelist files
! ----------------------------------------------

! 1.0 namelist nattrc :

      nwritetrc = 10
      lrsttr    =.FALSE.
      nrsttr    = 0

      DO ji=1,jptra
        WRITE (ctrcnm(ji),'("TR_",I1)') ji
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

      rsc     = 1.
      rtrn    = 1.e-15
      ncor    = 1
      ndttrc  = 4
      crosster= .FALSE.
      lhdf    = .TRUE.
      ahtrb0  = 0.
      trcrat  = 1.
      ahtrc0  = aht0
      aeivtr0 = aeiv0

      REWIND(numnat)
      READ(numnat,natnum)

      IF(lwp) THEN
          WRITE(numout,*) ' '
          WRITE(numout,*) 'natnum'
          WRITE(numout,*) ' '
          WRITE(numout,*) ' rsc tuning coefficient                  = ',rsc
          WRITE(numout,*) ' rtrn truncation value                   = ',rtrn
          WRITE(numout,*) ' ncor number of corrective phase         = ',ncor
          WRITE(numout,*) ' ndttrc time step freq. for pass. trac.  = ',ndttrc
          WRITE(numout,*) ' crosster  computes or not crossterms    = ',crosster
          WRITE(numout,*) ' lhdf  calls or not trchdf               = ',lhdf
          WRITE(numout,*) ' background diffusivity for passive tr   = ',ahtrb0
          WRITE(numout,*) ' ratio betweeen passive and active tr diffusion coeff= ' ,trcrat
          WRITE(numout,*) ' horizontal eddy diffus. for passive tr  = ',ahtrc0
          WRITE(numout,*) ' eddy induced veloc. coef. for passive tr = ',aeivtr0
          WRITE(numout,*) ' '
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

! 1.4 namelist natdmp : passive tracers dumping

!      CALL trclsm

      END SUBROUTINE trclec
