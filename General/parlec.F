! $Header: /cvsroot/opatm-bfm/opa_model/OPA/parlec.F,v 1.3 2010-03-05 09:56:44 cvsogs01 Exp $

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
!
!   Output :
!   -------
!      common
!            /comdom/          : domain parameters
!                                and scale factors
!            /comcoh/          : horizontal coordinates and factors
!            /comask/          : masks, bathymetry
!            /comhdt/          : lateral diffusivity
!                                diffusion
!            /comeos		: ocean physical parameters
!            /comstp/          : time step
!            /comrst/          : restart parameter
!            /comctl/          : parameters for the control
!            /comwri/          : parameters for output
!            /comphy/          : parameters for physical process
!            /comtsk/          : multitasking
!            /comdyn/          : dynamic data selection
!


       USE myalloc
       IMPLICIT NONE

! local declarations
! ==================
      INTEGER ji,jj,jk
      INTEGER iost, ivide1, ivide2
      INTEGER ilu, iused(1,100)
      INTEGER ilglo, ildta, ibloc, ilseq

      CHARACTER*32 clname
      CHARACTER*21 clold, clfor, clseq, clnew, cldir, clunf


      NAMELIST/namrun/ lrstar, calendarType, nsptint

      NAMELIST/namdom/ nizoom, njzoom, nacc, atfp, rdt, rdtmin, rdtmax, rdth
      NAMELIST/namhdf/ aht0, ahtb0, aeiv0
      NAMELIST/namzdf/ avt0,ntrbbl,atrbbl
      NAMELIST/nameos/ neos, rau0, ralpha, rbeta, rcp



!----------------------------------------------------------------------
! statement functions
! ===================

#include "stafun.h"
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
      REWIND( numnam )
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




      CLOSE( numnam)
      CALL FLUSH_(numout)
      END SUBROUTINE parlec
