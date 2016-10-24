            SUBROUTINE parcst
!!!---------------------------------------------------------------------
!!!
!!!                       ROUTINE parcst
!!!                     ******************
!!!
!!!  Purpose :
!!!  --------
!!!     Print the model parameters-
!!!----------------------------------------------------------------------

       USE myalloc
        IMPLICIT NONE


      if (lwp) then
      WRITE(numout,*) ' '
      WRITE(numout,*) ' routine parcst'
      WRITE(numout,*) ' **************'
      WRITE(numout,*) ' '
      WRITE(numout,*) ' initialization of ogstm model'
      WRITE(numout,*) ' '

! 1. Parameters
! -------------

      WRITE(numout,*) ' '
      WRITE(numout,*) ' parameter file'
      WRITE(numout,*) ' **************'
      WRITE(numout,*) ' '

! 1.1 Dimensions of domain

      WRITE(numout,*) ' dimension of model'
      WRITE(numout,*) ' '
      WRITE(numout,*) ' jpni    : ', jpni
      WRITE(numout,*) ' jpnj    : ', jpnj
      WRITE(numout,*) ' jpnij   : ', jpnij
      WRITE(numout,*) ' jpreci  : ', jpreci
      WRITE(numout,*) ' jprecj  : ', jprecj
      WRITE(numout,*) ' '
      WRITE(numout,*) ' jpiglo  : ', jpiglo
      WRITE(numout,*) ' jpjglo  : ', jpjglo
      WRITE(numout,*) ' jpk     : ', jpk
      WRITE(numout,*) ' '
      WRITE(numout,*) ' jpi     : ', jpi
      WRITE(numout,*) ' jpj     : ', jpj
      WRITE(numout,*) ' '
      WRITE(numout,*) ' jpim1   : ', jpim1
      WRITE(numout,*) ' jpjm1   : ', jpjm1
      WRITE(numout,*) ' jpkm1   : ', jpkm1
      WRITE(numout,*) ' jpij    : ', jpij
      WRITE(numout,*) ' '

! 1.2 Option parameter : island, periodic,...

      WRITE(numout,*) ' island,',' lateral boundary conditions, s-coordinate'
      WRITE(numout,*) ' '
      WRITE(numout,*) ' jperio  : ', jperio
      WRITE(numout,*) ' '


      endif

      END SUBROUTINE parcst



      SUBROUTINE parctl 
!!!---------------------------------------------------------------------
!!!
!!!                       ROUTINE parctl
!!!                     ******************
!!!
!!!  Purpose :
!!!  --------
!!!     Control the cpp options for the run and if files are availables
!!!     Control also consistancy between options and namelist values
!!!
!!   Method :
!!   -------
!!      We use if/endif inside #if defined option-cpp
!!      c a u t i o n : FILE name must not exceed 21 characters

       USE myalloc
       IMPLICIT NONE
!! local declarations
!! ==================
      INTEGER :: ji, jj
      INTEGER :: iadv, istop, iwarn



      IF(lwp) THEN
          WRITE(numout,*) ' '
          WRITE(numout,*) ' '
          WRITE(numout,*) ' routine parctl'
          WRITE(numout,*) ' **************'
          WRITE(numout,*) ' '
          WRITE(numout,*) ' control of cpp options and files'
          WRITE(numout,*) ' '
      END IF
!! Initialization
      istop = 0
      iwarn = 0
9000  FORMAT( /,' ===>>>> : E R R O R',/,'          ===========',/ )

!! II. Domain
!! ... Intialization of sub domain index
      DO ji = 1,jpi
        mindi(ji) = ji+nizoom-1
      END DO
      DO jj = 1,jpj
        mindj(jj) = jj+njzoom-1
      END DO

      IF (lwp) THEN
       WRITE(numout,*) ' '
       WRITE(numout,*) '*** domain '
       WRITE(numout,*) '    global domain    : jpiglo = ', jpiglo,' jpjglo = ', jpjglo, ' jpk    = ', jpk
       WRITE(numout,*) '    local domain     : jpi    = ', jpi,' jpj    = ', jpj   , ' jpk    = ', jpk
       WRITE(numout,*) '    index in i coordinate '
       WRITE(numout,25) (mindi(ji),ji = 1,jpi)
       WRITE(numout,*) ' '
       WRITE(numout,*) '    index in j coordinate '
       WRITE(numout,25) (mindj(jj),jj = 1,jpj)
       WRITE(numout,*) ' '
      END IF
 25   FORMAT( (100(4x,19i4,/)) )
      print *,"jp",jpi,jpj
      IF ( mindi(jpi).GT.jpiglo .OR. mindj(jpj).GT.jpjglo ) THEN
          IF(lwp)WRITE(numout,9000)
          IF(lwp)WRITE(numout,*) ' subdomain greater than the initial'
          IF(lwp)WRITE(numout,*) ' one, check your dimensions'
          istop = istop + 1
      ENDIF


!!... coordinates & bathymetry
!! ----------------------------
!!   ... Vertical coordinate

      IF(lwp) THEN
          WRITE(numout,*)
          WRITE(numout,*)
          WRITE(numout,*) ' *** vertical coordinate option'
          WRITE(numout,*) '     z-coordinates (default)'
!! VI. Data and surface forcing dynamic and physic input file
          WRITE(numout,*)
          WRITE(numout,*)
          WRITE(numout,*) ' *** dynamic and physics file'
          WRITE(numout,*) ' *** we USE temp. and salinity fields'
          WRITE(numout,*)

          WRITE(numout,*)
          WRITE(numout,*)
          WRITE(numout,*) '     Algorithmic and numerical schemes'
          WRITE(numout,*) '     ---------------------------------'
          WRITE(numout,*)
          WRITE(numout,*)
          WRITE(numout,*) '     Ocean Physics'
          WRITE(numout,*) '     -------------'
          WRITE(numout,*)
          WRITE(numout,*)
          WRITE(numout,*)
          WRITE(numout,*) ' *** equation of state option'
      ENDIF
      IF ( neos.EQ.0 ) THEN
          IF(lwp)WRITE(numout,*) '     use of Jackett & McDougall (1994) equation of state and'
          IF(lwp)WRITE(numout,*) '            McDougall (1987) Brunt-Vaisala frequency'
        ELSEIF ( neos.EQ.1 ) THEN
          IF(lwp)WRITE(numout,*) '     use of linear eos rho(T) = rau0 * ( 1.028 - ralpha * T )'
        ELSEIF ( neos.EQ.2 ) THEN
          IF(lwp)WRITE(numout,*) '     use of linear eos rho(T,S) = rau0 * ( rbeta * S - ralpha * T )'
        ELSE
          IF(lwp)WRITE(numout,9000)
          IF(lwp)WRITE(numout,*) ' neos flag has a wrong value : ',neos
          istop = istop + 1
      ENDIF



      IF(lwp) THEN
          WRITE(numout,*)
          WRITE(numout,*)
          WRITE(numout,*) ' *** vertical diffusion option'
          WRITE(numout,*) '     use an implicit scheme'
          WRITE(numout,*)
          WRITE(numout,*)
          WRITE(numout,*) ' *** lateral diffusion option'
          WRITE(numout,*)
!   ... Space variation of eddy coefficients
#if defined key_trahdfcoef1d 
          WRITE(numout,*) '     tracer eddy coef. function of depth only'
#endif
!   ... Type of diffusive operator
#if defined key_trahdfbilap
          WRITE(numout,*) '     biharmonic tracer diffusion'
#  else
          WRITE(numout,*) '     harmonic tracer diffusion (default)'
#endif
      ENDIF
!   ... Order logical units by growing numbers


! ... W a r n i n g  and  e r r o r  control
! ------------------------------------------

      IF ( istop.GT.0 ) THEN
          IF(lwp)WRITE(numout,*)
          IF(lwp)WRITE(numout,*) istop,' E R R O R found : we stop'
          IF(lwp)WRITE(numout,*) '**************************'
          IF(lwp)WRITE(numout,*)
          STOP 'parctl'
      ENDIF




      istop = 0

! 8. Advection scheme option
! --------------------------
      iadv = 0
      IF(lwp) THEN
          WRITE(numout,*) ' *** Advection scheme'
          WRITE(numout,*)
#if defined key_trc_smolar
          WRITE(numout,*) ' 3D advection with Smolarkiewicz scheme'
          WRITE(numout,*) ' '
          iadv = iadv + 1
#endif
! 9. Lateral diffusion option
          WRITE(numout,*)' *** Lateral diffusion option for passive tracer'
          WRITE(numout,*)
!   ... Type of diffusive operator
#     if defined key_trc_hdfbilap
         WRITE(numout,*) '   biharmonic tracer diffusion'
#     else
#        if defined key_trc_hdflap
         WRITE(numout,*) '   harmonic tracer diffusion'
#        else
         WRITE(numout,*) '   passive tracer diffusion (default)'
         WRITE(numout,*) '   samethan active tracer diffusion'
#        endif
#     endif
!     10. tracer damping option
          WRITE(numout,*) ' *** Tracer damping option'
          WRITE(numout,*)
#     if defined key_trc_dmp
          WRITE(numout,*)'key_trc_dmp is defined'
#    else
          WRITE(numout,*) ' No tracer damping'
#    endif
          WRITE(numout,*) ' *** Source/Sink model option'
          WRITE(numout,*)
! 11. SMS model
#    if defined key_trc_bfm
          WRITE(numout,*) ' use bfm tracer model '
#    else
          WRITE (numout,*) ' No Source/Sink model '
#    endif
          WRITE(numout,*) ' '




!      E r r o r  control

      IF ( istop.GT.0 .or. iadv.NE.1 ) THEN
          IF (iadv.EQ.0) THEN
             WRITE(numout,*) 'No advection scheme is defined'
             WRITE(numout,*) '******************************'
             istop = istop + 1
          ELSE
             WRITE(numout,*) iadv, 'advection schemes are defined'
             WRITE(numout,*) '***********************************'
             istop = istop + 1
          ENDIF
          IF(lwp)WRITE(numout,*)
          IF(lwp)WRITE(numout,*) istop,' E R R O R found : we stop'
          IF(lwp)WRITE(numout,*) '**************************'
          IF(lwp)WRITE(numout,*)
          STOP 'trcctl'
      ENDIF

      ENDIF


      END SUBROUTINE PARCTL
