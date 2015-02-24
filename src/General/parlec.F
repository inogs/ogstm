      SUBROUTINE parlec
!---------------------------------------------------------------------
!
!                       ROUTINE parlec
!                     ******************
!
!  Purpose :
!  --------
!     Read and print options for the run (namelist)
!
!   Method :                   : no
!   -------
!
!   Input :
!   ------
!      the namelist file ( unit numnam ) :
!            &namrun           : parameters of the run
!            &namdom           : domain and time parameters
!            &namhdf           : horizontal diffusion
!            &namzdf           : vertical diffusion
!            &nameos           : ocean physical parameters
!            &namdyn           : offline dynamic selection
!            &natnum           : numerical schemes
!            &general_IO       : IO settings



       USE myalloc
       IMPLICIT NONE

! local declarations
! ==================
      INTEGER ji!,jj,jk
!      INTEGER iost, ivide1, ivide2
      INTEGER ilu, iused(1,100)
      INTEGER ilglo, ildta, ibloc, ilseq

      CHARACTER(LEN=32) clname
      CHARACTER(LEN=21) clold, clfor, clseq, clnew, cldir, clunf


      NAMELIST/namrun/ lrstar, calendarType, nsptint

      NAMELIST/namdom/ nizoom, njzoom, nacc, atfp, rdt, rdtmin, rdtmax, rdth
      NAMELIST/namhdf/ aht0, ahtb0, aeiv0
      NAMELIST/namzdf/ avt0,ntrbbl,atrbbl
      NAMELIST/nameos/ neos, rau0, ralpha, rbeta, rcp
      namelist /natnum/ rsc,rtrn,ncor,ndttrc,crosster,lhdf,ahtrb0,trcrat,ahtrc0,aeivtr0,vsed,photop,atlantic_bfm,bottom_flux
!      namelist /natnum/ rsc,rtrn,ndttrc,crosster,lhdf,vsed,photop,atlantic_bfm,bottom_flux
!      namelist /dovesiete/  ncor,ahtrb0,trcrat,ahtrc0,aeivtr0
      NAMELIST/General_IO/   nwritetrc,lrsttr, nrsttr,freq_ave_phys,save_bkp_group2, isCheckLOG



! ... Open specifier
      clold = 'OLD'
      clnew = 'NEW'
      clfor = 'FORMATTED'
      clunf = 'UNFORMATTED'
      clseq = 'SEQUENTIAL'
      cldir = 'DIRECT'


      ilseq = 1 ! ... Sequential value

! ... computation of the record length for direct access file
!     this length depend of 512 for the t3d machine
      ibloc = 512
      ildta = ibloc*((jpidta*jpjdta)/ibloc+1)*jpbyt
      ilglo = ibloc*((jpiglo*jpjglo)/ibloc+1)*jpbyt
      ilu = 0
      DO ji = 1, 100
        iused(1,ji) = 0
      END DO

      IF(lwp) THEN
          WRITE(numout,*) ' '
          WRITE(numout,*) ' '
          WRITE(numout,*) ' routine parlec'
          WRITE(numout,*) ' **************'
          WRITE(numout,*) ' '
          WRITE(numout,*) ' initialization of run'
          WRITE(numout,*) ' '
          WRITE(numout,*) ' '
          WRITE(numout,*) ' namelist'
          WRITE(numout,*) ' ********'
          WRITE(numout,*) ' '
      ENDIF

      numnam = 4
      clname = 'namelist'
      CALL ctlopn( numnam, clname, clold, clfor, clseq,ilseq, ilu, iused, numout, lwp, 1)


! 1. Read namlist files
! ---------------------

! ... Namelist namdom : space/time domain (bathymetry, mesh, timestep)

      nizoom =  1
      njzoom =  1

      nacc   = 0
      atfp   = 0.1
      rdt    = 3600.
      rdtmin = 3600.
      rdtmax = 3600.
      rdth   = 800.

      REWIND( numnam )
      READ  ( numnam, namdom )

      IF(lwp) THEN
      WRITE(numout,*) ' '
      WRITE(numout,*) 'namdom'
      WRITE(numout,*) ' '
      WRITE(numout,*) ' index i domain left bottom   nizoom   = ',nizoom
      WRITE(numout,*) ' index j domain left bottom   njzoom   = ',njzoom
      WRITE(numout,*) ' '
      WRITE(numout,*) ' acceleration of converge       nacc   = ',nacc
      WRITE(numout,*) ' asselin time filter parameter  atfp   = ',atfp
      WRITE(numout,*) ' time step                      rdt    = ',rdt
      WRITE(numout,*) ' minimum time step on tracers   rdtmin = ',rdtmin
      WRITE(numout,*) ' maximum time step on tracers   rdtmax = ',rdtmax
      WRITE(numout,*) ' depth variation tracer step    rdth   = ',rdth
      WRITE(numout,*) ' '

      ENDIF

! ... Namelist namhdf : horizontal diffusion

      aht0   = 2000.
      ahtb0  = 0.
      aeiv0  = 2000.

      REWIND( numnam )
      READ  ( numnam, namhdf )

      IF(lwp) THEN
      WRITE(numout,*) 'namhdf'
      WRITE(numout,*) ' '
      WRITE(numout,*) ' horizontal eddy diffusivity    aht0   = ',aht0
      WRITE(numout,*) ' background diffusivity         ahtb0  = ',ahtb0
      WRITE(numout,*) ' eddy induced velocity coef.    aeiv0  = ',aeiv0
      WRITE(numout,*) ' '
      ENDIF

! ... Namelist namzdf : vertical diffusion

      avt0  = 1.e-5
      ntrbbl = 0
      atrbbl = 10000.

      REWIND( numnam )
      READ  ( numnam, namzdf )

      IF(lwp) THEN
      WRITE(numout,*) 'namzdf'
      WRITE(numout,*) ' '
      WRITE(numout,*) ' vertical eddy diffusivity      avt0   = ',avt0
      WRITE(numout,*) ' '
      WRITE(numout,*) ' bottom boundary layer ntrbbl   = ',ntrbbl
      WRITE(numout,*) ' '
      WRITE(numout,*) ' bottom boundary layer diff. atrbbl   = ',atrbbl
      WRITE(numout,*) ' '
      ENDIF


! ... namelist nameos : ocean physical parameters

      neos   =      0
      rau0   =  1020.
      ralpha =  2.e-4
      rbeta  =  0.029
      rcp    =  4.e+3

      REWIND( numnam )
      READ  ( numnam, nameos )

      ro0cpr = 1. / ( rau0 * rcp )

      IF(lwp) THEN
      WRITE(numout,*) 'nameos'
      WRITE(numout,*) ' '
      WRITE(numout,*) ' flag for eq. of state and N^2  neos   = ',neos
      WRITE(numout,*) ' volumic mass of reference      rau0   = ',rau0
      WRITE(numout,*) ' thermal exp. coef. (linear)    ralpha = ',ralpha
      WRITE(numout,*) ' saline exp. coef. (linear)     rbeta  = ',rbeta
      WRITE(numout,*) ' ocean specific heat            rcp    = ',rcp
      WRITE(numout,*) '          1. / ( rau0 * rcp ) = ro0cpr = ',ro0cpr
      WRITE(numout,*) ' '
      ENDIF

! ... Namelist namrun : parameters of the run


      lrstar = .FALSE.



      REWIND( numnam )
      READ  ( numnam, namrun )

      IF(lwp) THEN
      WRITE(numout,*) 'namrun'
      WRITE(numout,*) ' '
      WRITE(numout,*) ' restart logical                lrstar = ',lrstar
      WRITE(numout,*) ' leap year calendar (0/1)  calendarType = ',calendarType
      WRITE(numout,*) ' type of interpolation                     nsptint = ', nsptint

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

      REWIND(numnam)
      READ(numnam,natnum)

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



! ************* namelist GENERAL_IO  *****************

      nwritetrc = 10
      lrsttr    =.FALSE.
      nrsttr    = 0

      REWIND( numnam )
      READ  ( numnam, General_IO )



      if (lwp) then
          if (freq_ave_phys.eq.2) then
               WRITE(numout,*) 'Forcings phys will be dumped following 2.aveTimes file'
               else
               if (freq_ave_phys.eq.1) then
                 WRITE(numout,*) 'Forcings phys will be dumped following 1.aveTimes file'
               else
                 WRITE(numout,*) 'Forcings phys will not be dumped'
               endif
          endif
      endif



      CLOSE( numnam)

      END SUBROUTINE parlec
