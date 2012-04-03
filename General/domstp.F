
      SUBROUTINE domstp
!---------------------------------------------------------------------
!
!                       ROUTINE domstp
!                     ******************
!
!  Purpose :
!  --------
!     Intialize time step for the run
!
!   Method :
!   -------
!	Initialization of a coef. use in the Asselin time filter
!		atfp1 = 1 - 2 * atfp
!	where atfp is the Asselin time filter parameter read in namelist
!
!	Model time step:
!         nacc = 0 : synchronous time intergration.
!	There is one time step only, defined by: rdt, rdttra(k)=rdt
!         If nacc = 1 accelerating the convergence. There is 2
!      different time steps for dynamics and tracers:
!             		rdt      : dynamical part
!             		rdttra(k): temperature and salinity
!      The tracer time step is a function of vertical level. the model
!      reference time step ( i.e. for wind stress, surface heat and
!      salt fluxes) is the surface tracer time step is rdttra(1).
!         N.B. depth dependent acceleration of convergence is not im-
!      plemented for s-coordinate.
!
!   Input :
!   ------
!      common
!            /comstp/ nacc     : flag for the convergence acceleration
!                     rdt      : time step for dynamics
!                     rdtmin   : minimum time step on tracers
!                     rdtmax   : maximum time step on tracers
!                     rdth     : depth scale for vertically dependent
!                                time step for tracers
!
!   Output :
!   -------
!      common
!            /comstp/ rdttra   : vertical profile of tracer time step
!		       atfp1    : = 1 - 2*atfp
!
!   References :
!   -----------
!	Bryan, K., 1984, J. Phys. Oceanogr., 14, 666-673.
!
!   Modifications :
!   --------------
!      original :  90-10 (O. Marti)
!      additions : 96-01 (G. Madec) terrain following coordinates
!----------------------------------------------------------------------

       USE myalloc
       USE stafun
       IMPLICIT NONE

!    local declarations
      INTEGER jk

! statement functions
! ===================



! 0. Asselin Time filter
! ----------------------

      atfp1 = 1. - 2. * atfp


      IF ( nacc.EQ.0 ) THEN

! 1. Synchronous integration
! --------------------------

          IF(lwp)WRITE(numout,*)
          IF(lwp)WRITE(numout,*)' domstp: synchronous integration'
          IF(lwp)WRITE(numout,*)' ******* -----------------------'
          IF(lwp)WRITE(numout,*)
          IF(lwp)WRITE(numout,*)'     dynamics and tracer time step = ',rdt/3600., ' hours'

          DO jk = 1, jpk
            rdttra(jk) = rdt
          END DO  


        ELSEIF ( nacc.EQ.1 ) THEN

! 2. Accelerating the convergence
! -------------------------------

          IF(lwp)WRITE(numout,*)
          IF(lwp)WRITE(numout,*)' domstp: accelerating the convergence'
          IF(lwp)WRITE(numout,*)' ******* ----------------------------'
          IF(lwp)WRITE(numout,*)
          IF(lwp)WRITE(numout,*)'         dynamical time step : ',rdt/3600.,' hours'
#if defined key_s_coord
          IF ( rdtmin.NE.rdtmax ) THEN
              IF(lwp)WRITE(numout,*)'         e r r o r'
              IF(lwp)WRITE(numout,*)'         ========='
              IF(lwp)WRITE(numout,*)' depth dependent acceleration of convergence not implemented in s-coordinates'
              STOP 'domstp'
          ENDIF
#endif
          IF(lwp)WRITE(numout,*)'         tracers   time step :','  dt (hours)  level'

          DO jk = 1, jpk
            IF ( fsdept(1,1,jk).LE.rdth ) rdttra(jk) = rdtmin
            IF ( fsdept(1,1,jk).GT.rdth )
     $          rdttra(jk) = rdtmin
     $               + ( rdtmax - rdtmin )
     $               * ( exp( ( fsdept(1,1,jk ) - rdth ) / rdth ) - 1. )
     $               / ( exp( ( fsdept(1,1,jpk) - rdth ) / rdth ) - 1. )
            IF(lwp)WRITE(numout,9200) rdttra(jk)/3600., jk
          END DO  
 9200     FORMAT(36x,f5.2,'     ',i3)


      ELSE

! 3. E R R O R
! ------------

          IF(lwp)WRITE(numout,*)
          IF(lwp)WRITE(numout,*) ' nacc value e r r o r, nacc= ',nacc
          IF(lwp)WRITE(numout,*) ' we stop'
          STOP 'domstp'

      ENDIF


      END SUBROUTINE domstp
