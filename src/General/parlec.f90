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
!            &namhdf           : horizontal diffusion
!            &namzdf           : vertical diffusion
!            &nameos           : ocean physical parameters
!            &natnum           : numerical schemes
!            &general_IO       : IO settings



       USE myalloc
       ! epascolo USE myalloc_mpp
       IMPLICIT NONE

! local declarations
! ==================

      NAMELIST/namrun/ calendarType, nsptint

      NAMELIST/namhdf/ aht0
      NAMELIST/nameos/ neos, rau0, ralpha, rbeta
      namelist /natnum/ rdt,rsc,rtrn,ncor,ndttrc,lhdf,lrivers,lbfm, latmosph, ahtrb0,trcrat,ahtrc0,vsed,photop,atlantic_bfm,bottom_flux
      NAMELIST/General_IO/   nwritetrc, freq_ave_phys,save_bkp_group2, isCheckLOG

      NAMELIST/Domain_Characteristic/  jperio
      NAMELIST/Number_Fluxes/ jpflx, jpwind, jpemp,jpkef, jpice, jpqsr



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
      OPEN(unit=numnam, file='namelist.init', status='OLD')!, 'FORMATTED', 'SEQUENTIAL')

! 1. Read namlist files
! ---------------------



! ... Namelist namhdf : horizontal diffusion

      aht0   = 2000.


      REWIND( numnam )
      READ  ( numnam, namhdf )

      IF(lwp) THEN
      WRITE(numout,*) 'namhdf'
      WRITE(numout,*) ' '
      WRITE(numout,*) ' horizontal eddy diffusivity    aht0   = ',aht0
      WRITE(numout,*) ' '
      ENDIF


! ... namelist nameos : ocean physical parameters

      neos   =      0
      rau0   =  1020.
      ralpha =  2.e-4
      rbeta  =  0.029

      REWIND( numnam )
      READ  ( numnam, nameos )


      IF(lwp) THEN
      WRITE(numout,*) 'nameos'
      WRITE(numout,*) ' '
      WRITE(numout,*) ' flag for eq. of state and N^2  neos   = ',neos
      WRITE(numout,*) ' volumic mass of reference      rau0   = ',rau0
      WRITE(numout,*) ' thermal exp. coef. (linear)    ralpha = ',ralpha
      WRITE(numout,*) ' saline exp. coef. (linear)     rbeta  = ',rbeta
      WRITE(numout,*) ' '
      ENDIF

! ... Namelist namrun : parameters of the run




      REWIND( numnam )
      READ  ( numnam, namrun )

      IF(lwp) THEN
      WRITE(numout,*) 'namrun'
      WRITE(numout,*) ' '
      WRITE(numout,*) ' leap year calendar (0/1)  calendarType = ',calendarType
      WRITE(numout,*) ' type of interpolation                     nsptint = ', nsptint

      ENDIF



! 1.1 namelist natnum :
      rdt         = 3600.
      rsc         = 1.
      rtrn        = 1.e-15
      ncor        = 1
      ndttrc      = 4.0
      lhdf        = .TRUE.
      ahtrb0      = 0.
      trcrat      = 1.
      ahtrc0      = aht0
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
          WRITE(numout,*) ' time step                      rdt                = ',rdt
          WRITE(numout,*) ' rsc tuning coefficient                            = ' , rsc
          WRITE(numout,*) ' rtrn truncation value                             = ' , rtrn
          WRITE(numout,*) ' ncor number of corrective phase                   = ', ncor
          WRITE(numout,*) ' ndttrc time step freq. for pass. trac.            = ', ndttrc
          WRITE(numout,*) ' lhdf  calls or not trchdf                         = ', lhdf
          WRITE(numout,*) ' activation of rivers                              = ', lrivers
          WRITE(numout,*) ' activation atmospheric deposition                 = ', latmosph
          WRITE(numout,*) ' activation of bfm                                 = ', lbfm
          WRITE(numout,*) ' background diffusivity for passive tr             = ', ahtrb0
          WRITE(numout,*) ' ratio betweeen passive and active tr diffusion coeff= ', trcrat
          WRITE(numout,*) ' horizontal eddy diffus. for passive tr            = ', ahtrc0
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




! ... namelist.init: Domain_Characteristic

      REWIND( numnam )
      READ  ( numnam, Domain_Characteristic )


      IF(lwp) THEN
      WRITE(numout,*) 'Domain_Characteristic'
      WRITE(numout,*) ' '
      WRITE(numout,*) ' lateral cond. type for the global domain (4, 3, 2, 1 or 0)',jperio
      ENDIF

! ..  Namelist :  Number of fluxes


      REWIND( numnam )
      READ  ( numnam, Number_Fluxes )

      IF(lwp) THEN
      WRITE(numout,*) 'Number of fluxes'
      WRITE(numout,*) ' '
      WRITE(numout,*) ' jpflx   : number of fluxes',jpflx
      WRITE(numout,*) ' jpemp   :  E - P in mm/day',jpemp
      WRITE(numout,*) ' jpqsr   :  solar radiation',jpqsr
      WRITE(numout,*) ' jpkef   :  Extinction Coefficient',jpkef
      WRITE(numout,*) ' '
      ENDIF




      CLOSE( numnam)

      END SUBROUTINE parlec
