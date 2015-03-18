      SUBROUTINE parctl 
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE parctl
CCC                     ******************
CCC
CCC  Purpose :
CCC  --------
CCC     Control the cpp options for the run and if files are availables
CCC     Control also consistancy between options and namelist values
CCC
CC   Method :
CC   -------
CC      We use if/endif inside #if defined option-cpp
CC      c a u t i o n : FILE name must not exceed 21 characters

       USE myalloc
       IMPLICIT NONE
CC local declarations
CC ==================
      INTEGER ji, jj
      INTEGER istop, iwarn



      IF(lwp) THEN
          WRITE(numout,*) ' '
          WRITE(numout,*) ' '
          WRITE(numout,*) ' routine parctl'
          WRITE(numout,*) ' **************'
          WRITE(numout,*) ' '
          WRITE(numout,*) ' control of cpp options and files'
          WRITE(numout,*) ' '
      END IF
C Initialization
      istop = 0
      iwarn = 0
9000  FORMAT( /,' ===>>>> : E R R O R',/,'          ===========',/ )

C II. Domain
C ... Intialization of sub domain index
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

      IF ( mindi(jpi).GT.jpiglo .OR. mindj(jpj).GT.jpjglo ) THEN
          IF(lwp)WRITE(numout,9000)
          IF(lwp)WRITE(numout,*) ' subdomain greater than the initial'
          IF(lwp)WRITE(numout,*) ' one, check your dimensions'
          istop = istop + 1
      ENDIF


C ... coordinates & bathymetry
C ----------------------------
C   ... Vertical coordinate

      IF(lwp) THEN
          WRITE(numout,*)
          WRITE(numout,*)
          WRITE(numout,*) ' *** vertical coordinate option'
          WRITE(numout,*) '     z-coordinates (default)'
C VI. Data and surface forcing dynamic and physic input file
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
C   ... Space variation of eddy coefficients
#if defined key_trahdfcoef1d 
          WRITE(numout,*) '     tracer eddy coef. function of depth only'
#endif
C   ... Type of diffusive operator
#if defined key_trahdfbilap
          WRITE(numout,*) '     biharmonic tracer diffusion'
#  else
          WRITE(numout,*) '     harmonic tracer diffusion (default)'
#endif
      ENDIF
C   ... Order logical units by growing numbers


C ... W a r n i n g  and  e r r o r  control
C ------------------------------------------

      IF ( istop.GT.0 ) THEN
          IF(lwp)WRITE(numout,*)
          IF(lwp)WRITE(numout,*) istop,' E R R O R found : we stop'
          IF(lwp)WRITE(numout,*) '**************************'
          IF(lwp)WRITE(numout,*)
          STOP 'parctl'
      ENDIF

      END SUBROUTINE PARCTL
