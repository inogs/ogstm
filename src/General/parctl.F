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
CC      -------------
CC   External :
CC   ---------
CC      orders                   : sorting routine

       USE myalloc
       IMPLICIT NONE
CC local declarations
CC ==================
      INTEGER ji, jj,  ilu, ioptio
      INTEGER imulti,  imode, istop, iwarn
      INTEGER iused(1,100), iwork(257), indx(100)
      INTEGER ibloc, ilseq

      CHARACTER(LEN=21) clold, clfor, clseq, clnew, cldir, clunf, clunk

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
 9000 FORMAT( /,' ===>>>> : E R R O R',/,'          ===========',/ )
 9100 FORMAT( /,' ===>>>> : W A R N I N G',/,'          ===============',/ )

C I. Files and logical units
      ilu = 0
      DO ji = 1, 100
        iused(1,ji) = 0
      END DO
C ... Open specifier
      clold = 'OLD'
      clnew = 'NEW'
      clunk = 'UNKNOWN'
      clfor = 'FORMATTED'
      clunf = 'UNFORMATTED'
      clseq = 'SEQUENTIAL'
      cldir = 'DIRECT'

      ilseq = 1 ! Sequential value
C     computation of the record length for direct access file
C     this length depend of 512 for the t3d machine
      ibloc = 512

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

      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,*) ' *** vertical coordinate option'
      ioptio = 0
#ifdef key_s_coord
      IF(lwp)WRITE(numout,*) '     s-coordinates'
      ioptio = ioptio+1
#endif
      IF ( ioptio.eq.0) THEN
          IF(lwp)WRITE(numout,*) '     z-coordinates (default)'
      ENDIF

C VI. Data and surface forcing
C ... dynamic and physic input file

      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,*) ' *** dynamic and physics file'


#if defined key_off_tra
      IF(lwp) THEN
          WRITE(numout,*) ' *** we USE temp. and salinity fields'
          WRITE(numout,*)
      END IF
#endif


      IF(lwp) THEN
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
      ENDIF
C   ... Space variation of eddy coefficients
      ioptio = 0
#if defined key_trahdfcoef3d
      IF(lwp)WRITE(numout,*) '     tracer eddy coef. function of latitude, longitude and depth'
      ioptio = ioptio+1
#endif
#if defined key_trahdfcoef2d 
      IF(lwp)WRITE(numout,*) '     tracer eddy coef. function of latitude and longitude'
      ioptio = ioptio+1 
#endif
#if defined key_trahdfcoef1d 
      IF(lwp)WRITE(numout,*) '     tracer eddy coef. function of depth only'
      ioptio = ioptio+1 
#endif
      IF ( ioptio.EQ.0 ) THEN
          IF(lwp)WRITE(numout,*) '     constant tracer eddy coef. (default option)'
        ELSEIF ( ioptio.GT.1 ) THEN
          IF(lwp) THEN
          WRITE(numout,9000)
          WRITE(numout,*) ' use only one of the following keys: key_trahdfcoef3d, key_trahdfcoef2d, key_trahdfcoef1d'
          END IF
          istop = istop + 1
      ENDIF

#if defined key_s_coord
#    if defined key_trahdfcoef1d
      IF(lwp)WRITE(numout,9000)
      IF(lwp)WRITE(numout,*) ' key_trahdfcoef1d cannot be used with key_s_coord'
      istop = istop + 1
#    endif
#endif

C   ... Type of diffusive operator

#if defined key_trahdfbilap
      IF(lwp)WRITE(numout,*) '     biharmonic tracer diffusion'
#  else
      IF(lwp)WRITE(numout,*) '     harmonic tracer diffusion (default)'
#endif

C   ... direction of diffusion

      ioptio = 0
#if defined key_trahdfiso
      IF(lwp)WRITE(numout,*) '     isopycnal tracer diffusion'
      ioptio = ioptio+1
#endif
#if defined key_trahdfgeop
      IF(lwp)WRITE(numout,*) '     geopotential tracer diffusion'
#  if ! defined key_s_coord
      IF(lwp)WRITE(numout,9000)
      IF(lwp)WRITE(numout,*) ' key_trahdfgeop require key_s_coord'
      istop = istop + 1
#  endif
      ioptio = ioptio+1
#endif
      IF ( ioptio.EQ.0 ) THEN
          IF(lwp)WRITE(numout,*) '     iso-level tracer diff. (default)'
          IF(lwp)WRITE(numout,*)
        ELSEIF ( ioptio.GT.1 ) THEN
          IF(lwp)WRITE(numout,9000)
          IF(lwp)WRITE(numout,*) ' key_trahdfiso and key_trahdfgeop '
          IF(lwp)WRITE(numout,*) ' are incompatibles. '
          istop = istop + 1
      ENDIF
#if defined key_trahdfeiv
      IF(lwp)WRITE(numout,*) '     eddy induced velocity on tracers'
#  if ! defined key_trahdfiso
      IF(lwp)WRITE(numout,9000)
      IF(lwp)WRITE(numout,*) ' the eddy induced velocity on tracers'
      IF(lwp)WRITE(numout,*) ' requires isopycnal diffusion, defines'
      IF(lwp)WRITE(numout,*) ' key_trahdfiso. '
      istop = istop + 1
#  endif
#  if defined key_trahdfbilap
      IF(lwp)WRITE(numout,9000)
      IF(lwp)WRITE(numout,*) ' the eddy induced velocity on tracers'
      IF(lwp)WRITE(numout,*) ' does not work with biharmonic operator'
      IF(lwp)WRITE(numout,*) ' suppress key_trahdfbilap. '
      istop = istop + 1
#  endif
#endif



C ... FINAL CHECKS
C ... Logical units in use

      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,*) ' *** units opened for the run'
      IF(lwp)WRITE(numout,*)
C   ... Output file ( open in the ogstm routine )
      ilu = ilu+1
      iused(1,ilu) = numout

C   ... namelist file ( open in the parlec routine )
      ilu = ilu+1
      iused(1,ilu) = numnam

C   ... Order logical units by growing numbers
      imode = 0
      CALL orders(imode,iwork,iused,indx,ilu,1,8,1)
      IF ( imode.NE.0 ) THEN
          IF(lwp)WRITE(numout,9000)
          IF(lwp)WRITE(numout,*)' problem in sorting logical units'
          IF(lwp)WRITE(numout,*)' return code is: ',imode
          istop = istop + 1
      ENDIF
C   ... Print logical units in use
      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,*) '     Logical units used by the code :'
      IF(lwp)WRITE(numout,*) '     --------------------------------'
      IF(lwp)WRITE(numout,*) ( iused(1,indx(ji)), ji = 1, ilu )
      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,*) ' if you need a new unit, use any other number between 10 and 100'
C ... W a r n i n g  and  e r r o r  control
C ------------------------------------------
      IF ( iwarn.GT.0 ) THEN
          IF(lwp)WRITE(numout,*)
          IF(lwp)WRITE(numout,*) iwarn,' w a r n i n g found'
          IF(lwp)WRITE(numout,*) '**************************'
          IF(lwp)WRITE(numout,*)
      ENDIF
      IF ( istop.GT.0 ) THEN
          IF(lwp)WRITE(numout,*)
          IF(lwp)WRITE(numout,*) istop,' E R R O R found : we stop'
          IF(lwp)WRITE(numout,*) '**************************'
          IF(lwp)WRITE(numout,*)
          STOP 'parctl'
      ENDIF

      END SUBROUTINE PARCTL
