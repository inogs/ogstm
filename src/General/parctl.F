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
CC
CC   Input :                      no
CC   ------
CC
CC   Output :                     no
CC   -------
CC
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
      INTEGER ildta, ilglo, ibloc, ilseq, ilindta

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
C
      istop = 0
      iwarn = 0
 9000 FORMAT( /,' ===>>>> : E R R O R',/,'          ===========',/ )
 9100 FORMAT( /,' ===>>>> : W A R N I N G',/,'          ===============',/ )


#if defined restart || defined zdftke || defined bsfpcg || defined islands || defined qsr || defined saldta || defined temdta || defined tau || defined flx  || defined eoslinear || defined eoseckart || defined eosunesco || defined nobsf || defined bsfpcg || defined bsfsor || defined zdfconstant || defined zdfrichardson || defined temdmp || defined saldmp || defined dyndmp || defined zdfexplicit || defined npc || defined freeslip || defined dynhdfbilap || defined trahdfbilap || defined vorenergy || defined vorcombined || defined coupled
      IF(lwp)WRITE(numout,9000)
      IF(lwp)WRITE(numout,*) '  cpp keys have been changed '
      IF(lwp)WRITE(numout,*) ' All keynames begin with the prefix key_'
      IF(lwp)WRITE(numout,*) ' see the reference manual'
      istop = istop + 1
#endif

C
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C
C
C     --------------------------------------
C
#if defined key_trahdfisog
      IF(lwp)WRITE(numout,9000)
      IF(lwp)WRITE(numout,*) '  release: key_trahdfisog has been '
      IF(lwp)WRITE(numout,*) ' changed into key_trahdfeiv and requires '
      IF(lwp)WRITE(numout,*) ' key_trahdfiso to be defined.'
      istop = istop + 1
#endif
C
C
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C
C I. Files and logical units
C
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C
C
      ilu = 0
      DO ji = 1, 100
        iused(1,ji) = 0
      END DO
C
C ... Open specifier
C ------------------
C
      clold = 'OLD'
      clnew = 'NEW'
      clunk = 'UNKNOWN'
      clfor = 'FORMATTED'
      clunf = 'UNFORMATTED'
      clseq = 'SEQUENTIAL'
      cldir = 'DIRECT'
C
C     Sequential value
C
      ilseq = 1
C
C     computation of the record length for direct access file
C     this length depend of 512 for the t3d machine
C
      ibloc = 512
      ildta = ibloc*((jpidta*jpjdta)/ibloc+1)*jpbyt
      ilglo = ibloc*((jpiglo*jpjglo)/ibloc+1)*jpbyt
      ilindta = ibloc*((jpidta*jpjdta)/ibloc+1)*jpbytda
C
C ... Logical unit initialization for files
C
C     numout : output file             unit (always used)
C     numnam : namelist                unit (always used)
C     nummsh : coordinates/meshes      unit 
C     numbat : bathymetry              unit
C     numhba : bathymetry (meters)     unit (#ifdef "key_s_coord")
C     numdyn : dynamic and physic file unit
C     numtke : tke file                unit
C

C     numout = unit intialize in the main program ogstm < 10

C     numnam = unit intialize in parlec < 10
      nummsh = 13
C      numbat = 11 PAS DE BATHY.
      numdyn = 63
      numtke = 64
      numslp = 65
      numvgm = 66
C      numhba = 77

C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C
C II. Domain
C
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C
C
C ... Intialization of sub domain index
C -------------------------------------
C
      DO ji = 1,jpi
        mindi(ji) = ji+nizoom-1
      END DO
      DO jj = 1,jpj
        mindj(jj) = jj+njzoom-1
      END DO

      IF (lwp) THEN
       WRITE(numout,*) ' '
       WRITE(numout,*) '*** domain '
       WRITE(numout,*) '    data input domain: jpidta = ', jpidta,' jpjdta = ', jpjdta, ' jpkdta = ', jpkdta
       WRITE(numout,*) '    global domain    : jpiglo = ', jpiglo,' jpjglo = ', jpjglo, ' jpk    = ', jpk
       WRITE(numout,*) '    local domain     : jpi    = ', jpi,' jpj    = ', jpj   , ' jpk    = ', jpk
       WRITE(numout,*)
       WRITE(numout,*) '    south-west indices','    nizoom = ',nizoom,' njzoom = ',njzoom
       WRITE(numout,*) ' '
       WRITE(numout,*) '    index in i coordinate '
       WRITE(numout,25) (mindi(ji),ji = 1,jpi)
       WRITE(numout,*) ' '
       WRITE(numout,*) '    index in j coordinate '
       WRITE(numout,25) (mindj(jj),jj = 1,jpj)
       WRITE(numout,*) ' '
      END IF
 25   FORMAT( (100(4x,19i4,/)) )

      IF ( mindi(jpi).GT.jpidta .OR. mindj(jpj).GT.jpjdta ) THEN
          IF(lwp)WRITE(numout,9000)
          IF(lwp)WRITE(numout,*) ' subdomain greater than the initial'
          IF(lwp)WRITE(numout,*) ' one, check your dimensions'
          istop = istop + 1
      ENDIF


      IF(lwp)WRITE(numout,*) ' '
      IF(lwp)WRITE(numout,*) '**************************************** '
      IF(lwp)WRITE(numout,*) '* OPA8.1 OFFLINE                       * '
      IF(lwp)WRITE(numout,*) '*                                      * '
      IF(lwp)WRITE(numout,*) '*  global ocean simulation             * '
      IF(lwp)WRITE(numout,*) '* ------------------------             * '
      IF(lwp)WRITE(numout,*) '*                                      * '
      IF(lwp)WRITE(numout,*) '**************************************** '
      IF(lwp)WRITE(numout,*) '  '
      IF ( jpiglo.LT.jpidta .OR. jpjglo.LT.jpjdta ) THEN
          IF(lwp)WRITE(numout,*) ' '
          IF(lwp)WRITE(numout,*) '*********************************** '
          IF(lwp)WRITE(numout,*) '*                                 * '
          IF(lwp)WRITE(numout,*) '* use a subdomain zoom            * '
          IF(lwp)WRITE(numout,*) '*                                 * '
          IF(lwp)WRITE(numout,*) '*********************************** '
          IF(lwp)WRITE(numout,*) ' '
          IF(lwp)WRITE(numout,*) ' c a u t i o n'
          IF(lwp)WRITE(numout,*) ' ============='
          IF(lwp)WRITE(numout,*) 'don t use dist.coast file, compute it'
          IF(lwp)WRITE(numout,*) 'see island identification in dommba.F'
          IF(lwp)WRITE(numout,*) 'and diffusivity coef. in inimix.F'
          IF(lwp)WRITE(numout,*) 'and damping coef. in dtacof.F'
          IF(lwp)WRITE(numout,*) ' '
      ENDIF
C
C
C ... coordinates & bathymetry
C ----------------------------
C
C   ... Vertical coordinate
C
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


C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C
C VI. Data and surface forcing
C
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C
C
C ... dynamic and physic input file
C ---------------------------------
C
      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,*) ' *** dynamic and physics file'


#if defined key_off_tra
      IF(lwp) THEN
          WRITE(numout,*) ' *** we USE temp. and salinity fields'
          WRITE(numout,*)
      END IF
#endif
#if defined key_off_degrad && defined key_trahdfiso
C ... slope input file
C -------------------
C
      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,*) ' *** slope file'


      CALL ctlopn(numslp,'slp.dta',clold,clunf,clseq,ilindta,ilu,iused,numout,lwp,1)

#endif
#if defined key_off_degrad && defined key_trahdfeiv
C ... slope input file
C -------------------
C
      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,*) ' *** Gent et Mc William Speed  file'
C
C
      clname='vgm.dta'
      CALL ctlopn(numvgm,clname,clold,clunf,clseq,ilindta,ilu,iused,numout,lwp,1)

#endif

C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C
C Algorithmic & Schemes
C
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,*) '     Algorithmic and numerical schemes'
      IF(lwp)WRITE(numout,*) '     ---------------------------------'



#if defined key_mpp
      IF(lwp)WRITE(numout,*) '     use of mpp option (cray t3d+pvm): macrotasking is not used '
      IF ( jpni*jpnj.EQ.jpnij ) THEN
          IF(lwp)WRITE(numout,*) '     jpni x jpnj = jpnij : standard cutting out'
        ELSE
          IF ( jpni*jpnj.LT.jpnij ) THEN
              IF(lwp)WRITE(numout,9000)
              IF(lwp)WRITE(numout,*) ' jpnij > jpni x jpnj impossible'
              istop = istop + 1
          ENDIF
          IF(lwp)WRITE(numout,*) '     jpnij < jpni x jpnj : eliminate land processors algorithm'
      ENDIF

      ioptio = ioptio+1
#  else
      IF ( jpni.NE.1 .or. jpnj.NE.1 ) THEN
          IF(lwp)WRITE(numout,9000)
          IF(lwp)WRITE(numout,*) ' jpni x jpnj not equal to 1. The do-'
          IF(lwp)WRITE(numout,*) ' main is lay out to do mpp computer'
          istop = istop + 1
      ENDIF
#endif
      IF ( ioptio.EQ.0 ) THEN
          imulti = 1
          IF(lwp)WRITE(numout,*) '     macro-tasking option (default)'

          IF ( 1.GT.jpcpu ) THEN
              IF(lwp)WRITE(numout,9000)
              IF(lwp)WRITE(numout,*) ' ncpu=',1,' > ',' jpcpu=',jpcpu
              istop = istop + 1
          ENDIF
        ELSEIF ( ioptio.GT.1 ) THEN
          IF(lwp)WRITE(numout,9000)
          IF(lwp)WRITE(numout,*) ' key_monotasking and key_mpp options'
          IF(lwp)WRITE(numout,*) ' are incompatible'
          istop = istop + 1
      ENDIF
C
C
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C
C ... OCEAN PHYSICS
C
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C
      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,*) '     Ocean Physics'
      IF(lwp)WRITE(numout,*) '     -------------'
      IF(lwp)WRITE(numout,*)
C
C ... Equation of state
C ---------------------
C
      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,*) ' *** equation of state option'
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
C
C
C
C
C ... Vertical diffusion (implicit)
C ---------------------------------
C
      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,*) ' *** vertical diffusion option'
      IF(lwp)WRITE(numout,*) '     use an implicit scheme'
C
C
C ... Lateral Ocean Diffusion
C ---------------------------
C
      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,*) ' *** lateral diffusion option'
      IF(lwp)WRITE(numout,*)
C
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
C
#if defined key_s_coord
#    if defined key_trahdfcoef1d
      IF(lwp)WRITE(numout,9000)
      IF(lwp)WRITE(numout,*) ' key_trahdfcoef1d cannot be used with key_s_coord'
      istop = istop + 1
#    endif
#endif
C
C   ... Type of diffusive operator
C
#if defined key_trahdfbilap
      IF(lwp)WRITE(numout,*) '     biharmonic tracer diffusion'
#  else
      IF(lwp)WRITE(numout,*) '     harmonic tracer diffusion (default)'
#endif
C
C   ... direction of diffusion
C
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
C
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
C
C
C
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C
C ... FINAL CHECKS
C
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C
C
C
C ... Logical units in use
C ------------------------
C
      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,*) ' *** units opened for the run'
      IF(lwp)WRITE(numout,*)
C
C   ... Output file ( open in the ogstm routine )
      ilu = ilu+1
      iused(1,ilu) = numout

C
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
C
C   ... Print logical units in use
      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,*) '     Logical units used by the code :'
      IF(lwp)WRITE(numout,*) '     --------------------------------'
      IF(lwp)WRITE(numout,*) ( iused(1,indx(ji)), ji = 1, ilu )
      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,*) ' if you need a new unit, use any other number between 10 and 100'
C
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

      RETURN
      END
