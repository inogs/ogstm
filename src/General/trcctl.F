      SUBROUTINE trcctl
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE trcctl
CCC                     ******************
CCC
CCC  PURPOSE :
CCC  ---------
CCC     only for passive tracer
CCC     control the cpp options for the run and IF files are availables
CCC     control also consistancy between options and namelist values
CCC
CC   METHOD :
CC   -------
CC      we use IF/ENDIF inside #IF defined option-cpp
CC      c a u t i o n : FILE name must not exceed 21 characters
CC      -------------
CC
CC   INPUT :                      no
CC   -----
CC
CC   OUTPUT :                     no
CC   ------
CC
CC   WORKSPACE :
CC   ---------
CC      local
CC           clold,clnew,clfor,clunf,clseq,cldir,clname,
CC           ildta,ilglo
CC
CC   EXTERNAL :
CC   --------
CC
CC   MODIFICATIONS:
CC   --------------
CC      original :
CC                 04/00 (O. Aumont, M.A. Foujols) HAMOCC3 and P3ZD
CC      additions : 00/05 (A. Estublier) TVD Limiter Scheme
CC      additions : 00/06 (A. Estublier) MUSCL Scheme
CC      additions : 00/11 (MA Foujols, E Kestenare) Lateral diffusion option
CC      additions : 00/12  ( E Kestenare): improve controls of defined keys
CC 
       USE myalloc
       IMPLICIT NONE

CC local declarations
CC ==================

#if defined key_passivetrc
      CHARACTER(LEN=32) clname,clold,clfor,clseq,clnew,cldir,clunf,clunk
      INTEGER iused(1,100)
      INTEGER iadv,ildta,ilglo,ioptio,istop
      INTEGER jn


C 1. initialization
C -----------------
C
C 0. Parameter

      istop = 0

C 1. OPEN specifier

      clold='OLD'
      clnew='NEW'
      clunk='UNKNOWN'
      clfor='FORMATTED'
      clunf='UNFORMATTED'
      clseq='SEQUENTIAL'
      cldir='DIRECT'


C
C 3. LOGICAL UNIT initialization for specific files for passive tracer
C
C     nutrst : restart FILE  INPUT  UNIT (lrsttr=.TRUE.)
C     nutini(jptra) : UNIT for initial FILE for tracer
C
      nutrst = 83
      DO jn=1,jptra
        nutini(jn)= 0
      END DO 
C
C
C
C
C 8. Advection scheme option
C --------------------------
C
      iadv = 0
C
      IF(lwp) THEN
          WRITE(numout,*) ' *** Advection scheme'
          WRITE(numout,*)
      ENDIF

#if defined key_trc_arakawa
      IF(lwp) THEN
          WRITE(numout,*) ' horizontal and vertical advection '
          WRITE(numout,*) ' with Arakawa scheme'
          WRITE(numout,*) ' '
      ENDIF
      iadv = iadv + 1
#endif
#if defined key_trc_tvd
      IF(lwp) THEN
          WRITE(numout,*) ' horizontal and vertical advection '
          WRITE(numout,*) ' with TVD scheme'
          WRITE(numout,*) ' '
      ENDIF
      iadv = iadv + 1
#endif
#if defined key_trc_muscl
      IF(lwp) THEN
          WRITE(numout,*) ' horizontal and vertical advection '
          WRITE(numout,*) ' with MUSCL scheme'
          WRITE(numout,*) ' '
      ENDIF
      iadv = iadv + 1
#    if defined key_mpp
      IF (lwp) THEN
          WRITE (numout,*) ' ===>>>> : w a r n i n g '
          WRITE (numout,*) ' =======   ============= '
          WRITE (numout,*) ' STOP, this scheme is not validate'
          WRITE (numout,*) ' on mpp mode'
      ENDIF 
      istop = istop + 1
#    endif
#endif
#if defined key_trc_ppm
      IF(lwp) THEN
          WRITE(numout,*) ' horizontal and vertical advection '
          WRITE(numout,*) ' with PPM scheme'
          WRITE(numout,*) ' '
      ENDIF
      iadv = iadv + 1
#endif
#if defined key_trc_smolar
      IF(lwp) THEN
          WRITE(numout,*) ' 3D advection with Smolarkiewicz scheme'
          WRITE(numout,*) ' '
      ENDIF
      iadv = iadv + 1
#endif

C
C 9. Lateral diffusion option
C ---------------------------
C
      IF(lwp) THEN
          WRITE(numout,*)' *** Lateral diffusion option for passive tracer'
          WRITE(numout,*)
      ENDIF
C
C   ... Type of diffusive operator
C
#     if defined key_trc_hdfbilap
         IF(lwp)WRITE(numout,*) '   biharmonic tracer diffusion'
#     else
#        if defined key_trc_hdflap
         IF(lwp)WRITE(numout,*) '   harmonic tracer diffusion'
#        else
         IF(lwp)WRITE(numout,*) '   passive tracer diffusion (default)'
         IF(lwp)WRITE(numout,*) '   samethan active tracer diffusion'
#        endif
#     endif
C
C   ... direction of diffusion
C
      ioptio = 0
#if defined key_trc_hdfiso
      IF(lwp)WRITE(numout,*) '     isopycnal tracer diffusion'
      ioptio = ioptio+1
#endif
#if defined key_trc_hdfgeop
      IF(lwp)WRITE(numout,*) '     geopotential tracer diffusion'
#  if ! defined key_s_coord
      IF(lwp)WRITE(numout,9000)
      IF(lwp)WRITE(numout,*) ' key_trc_hdfgeop require key_s_coord'
      istop = istop + 1
#  endif
      ioptio = ioptio+1
#endif
      IF ( ioptio.EQ.0 ) THEN
          IF(lwp)WRITE(numout,*) '     iso-level tracer diff. (default)'
          IF(lwp)WRITE(numout,*)
        ELSEIF ( ioptio.GT.1 ) THEN
          IF(lwp)WRITE(numout,9000)
          IF(lwp)WRITE(numout,*) ' key_trc_hdfiso and key_trc_hdfgeop '
          IF(lwp)WRITE(numout,*) ' are incompatibles. '
          istop = istop + 1
      ENDIF
C
#if defined key_trc_hdfeiv
      IF(lwp)WRITE(numout,*) '     eddy induced velocity on tracers'
#  if ! defined key_trc_hdfiso
      IF(lwp)WRITE(numout,9000)
      IF(lwp)WRITE(numout,*) ' the eddy induced velocity on tracers'
      IF(lwp)WRITE(numout,*) ' requires isopycnal diffusion, defines'
      IF(lwp)WRITE(numout,*) ' key_trc_hdfiso. '
      istop = istop + 1
#  endif
#  if defined key_trc_hdfbilap
      IF(lwp)WRITE(numout,9000)
      IF(lwp)WRITE(numout,*) ' the eddy induced velocity on tracers'
      IF(lwp)WRITE(numout,*) ' does not work with biharmonic operator'
      IF(lwp)WRITE(numout,*) ' suppress key_trc_hdfbilap. '
      istop = istop + 1
#  endif
#endif

9000  FORMAT( /,' ===>>>> : E R R O R',/,'          ===========',/ )


C 10. tracer damping option
C ---------------------------------------------
C
      IF(lwp) THEN
          WRITE(numout,*) ' *** Tracer damping option'
          WRITE(numout,*)
      ENDIF
#    if defined key_trc_dmp
      IF(lwp) THEN 
          WRITE(numout,*)'key_trc_dmp is defined'
          WRITE(numout,*)'Check trcdta SUBROUTINE '
          WRITE(numout,*)' '
      ENDIF 
#    else
      IF (lwp) WRITE(numout,*) ' No tracer damping'
#    endif

C
C
C 11. SMS model
C ---------------------------------------------
C
      IF(lwp) THEN
          WRITE(numout,*) ' *** Source/Sink model option'
          WRITE(numout,*)
      ENDIF

#    if defined key_trc_npzd && defined key_trc_nnpzddom
      IF (lwp) THEN 
          WRITE (numout,*) ' ===>>>> : w a r n i n g '
          WRITE (numout,*) ' =======   ============= '
          WRITE (numout,*)' STOP, only one model can be specified '
      END IF 
      istop = istop + 1
#    endif
#    if defined key_trc_npzd && defined key_trc_hamocc3
      IF (lwp) THEN 
          WRITE (numout,*) ' ===>>>> : w a r n i n g '
          WRITE (numout,*) ' =======   ============= '
          WRITE (numout,*)' STOP, only one model can be specified '
      END IF 
      istop = istop + 1
#    endif
#    if defined key_trc_hamocc3 && defined key_trc_nnpzddom
      IF (lwp) THEN 
          WRITE (numout,*) ' ===>>>> : w a r n i n g '
          WRITE (numout,*) ' =======   ============= '
          WRITE (numout,*)' STOP, only one model can be specified '
      END IF 
      istop = istop + 1
#    endif

#    if defined key_trc_bfm
#     include "trcctl.bfm.h" 
#    else
      IF (lwp) THEN 
          WRITE (numout,*) ' No Source/Sink model '
          WRITE (numout,*) ' '
      END IF 
#    endif
C
C E r r o r  control
C ------------------
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

C
#else
C
C no passive tracers
C
#endif

      END SUBROUTINE TRCCTL
