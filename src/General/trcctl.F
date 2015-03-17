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
CC
CC   WORKSPACE :
CC   ---------
CC
CC   EXTERNAL :
CC   --------
       USE myalloc
       IMPLICIT NONE

CC local declarations
CC ==================

#if defined key_passivetrc
      INTEGER iadv,istop

      istop = 0

C 8. Advection scheme option
C --------------------------
      iadv = 0
      IF(lwp) THEN
          WRITE(numout,*) ' *** Advection scheme'
          WRITE(numout,*)
#if defined key_trc_smolar
          WRITE(numout,*) ' 3D advection with Smolarkiewicz scheme'
          WRITE(numout,*) ' '
          iadv = iadv + 1
#endif
C 9. Lateral diffusion option
          WRITE(numout,*)' *** Lateral diffusion option for passive tracer'
          WRITE(numout,*)
C   ... Type of diffusive operator
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
C     10. tracer damping option
          WRITE(numout,*) ' *** Tracer damping option'
          WRITE(numout,*)
#     if defined key_trc_dmp
          WRITE(numout,*)'key_trc_dmp is defined'
#    else
          WRITE(numout,*) ' No tracer damping'
#    endif
          WRITE(numout,*) ' *** Source/Sink model option'
          WRITE(numout,*)
C 11. SMS model
#    if defined key_trc_bfm
          WRITE(numout,*) ' use bfm tracer model '
#    else
          WRITE (numout,*) ' No Source/Sink model '
#    endif
          WRITE(numout,*) ' '

      ENDIF


C      E r r o r  control

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


#endif

      END SUBROUTINE TRCCTL
