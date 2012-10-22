CCC                       trclsm.sed.h
CCC                       *************
CCC
CCC  PURPOSE :
CCC  ---------
CCC     READs and PRINT options for BFM-like namelist
CCC
CC   METHOD :                   : no
CC   -------
CC
CC   INPUT :
CC   -----
CC            &natsed           : biological parameters 
CC
CC   OUTPUT :
CC   ------
CC      COMMON                  : (#ifdef "key_passivetrc")

CC            /cotsed/          : sink parameters

CC   WORKSPACE :                : no
CC   ---------
CC
CC   MODIFICATIONS:
CC   --------------
CC      original  : 99-10 (M.A. Foujols, M. Levy) passive tracer
CC      additions : 01-01 (E. Kestenare) add sediments
CC----------------------------------------------------------------------

CC----------------------------------------------------------------------
CC local declarations
CC ==================
CCC 17 11 2004 F79 already declared in trclsm.opt.h
CC    INTEGER ji
CC    INTEGER ilu,iused(1,100),ilseq
CC    CHARACTER(LEN=21) clold,clfor,clseq,clnew,cldir,clunf
CC    CHARACTER(LEN=32) clname

C 0. initializations
C ------------------
C
CCC 18 11 2004 F79 you can find it in trclsm.opt.h
CCC      namelist/natsed/ vsed,tmaxr,tminr
C
C
C
C
C initialize the number of LOGICAL UNIT used
C
      IF(lwp) THEN
          WRITE(numout,*) ' '
          WRITE(numout,*) ' ROUTINE trclec'
          WRITE(numout,*) ' **************'
          WRITE(numout,*) ' '
          WRITE(numout,*) ' namelist BFM for  model'
          WRITE(numout,*) ' ***********************'
          WRITE(numout,*) ' '
      ENDIF
C
C
C
C
C
C 1.4 namelist natsed : sinking  parameters
C
      vsed  = 0
      vsed_dia  = 0
      tmaxr   = 1./(     4.*rday)*0.
      tminr   = 1./(24.*30.*rday)*0.
C
C
      READ(numnat,natsed)
C
      IF(lwp) THEN
          WRITE(numout,*)
     $        ' detritus sedimentation speed          vsed =',86400*vsed
          WRITE(numout,*)
     $        ' maximum damping for d z or p         tmaxr =', tmaxr
          WRITE(numout,*)
     $        ' minimum damping for d z or p         tminr =', tminr
      ENDIF
