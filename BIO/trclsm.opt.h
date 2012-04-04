CCC                       trclsm.npzd.h
CCC                       *************
CCC
CCC  PURPOSE :
CCC  ---------
CCC     READs and PRINT options for NPZD namelist
CCC
CC   METHOD :                   : no
CC   -------
CC
CC   INPUT :
CC   -----
CC            &natbio           : biological parameters 
CC            &natopt           : optical parameters 
CC
CC   OUTPUT :
CC   ------
CC      COMMON                  : (#ifdef "key_passivetrc")

CC            /cotbio/          : biological parameters
CC            /cotopt/          : optical parameters

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

#if defined key_passivetrc && defined key_trc_npzd || defined key_trc_bfm 
      INTEGER ji
      INTEGER ilu,iused(1,100),ilseq
      CHARACTER(LEN=21) clold,clfor,clseq,clnew,cldir,clunf
      CHARACTER(LEN=32) clname

C 0. initializations
C ------------------
C
      namelist/natopt/xkg0,xkr0,xkgp,xkrp,xlg,xlr,rpig
      namelist/natsed/ vsed, vsed_dia, tmaxr, tminr
#if defined key_trc_diabio
      namelist/natdbi/ctrbio,ctrbil,ctrbiu,nwritebio
#endif

C
C
C OPEN specifier
C
      clold='old'
      clnew='new'
      clfor='formatted'
      clunf='unformatted'
      clseq='SEQUENTIAL'
      cldir='direct'
C
C SEQUENTIAL value
C
      ilseq=1
C
C initialize the number of LOGICAL UNIT used
C
      ilu=0

      IF(lwp) THEN
          WRITE(numout,*) ' '
          WRITE(numout,*) ' ROUTINE trclec'
          WRITE(numout,*) ' **************'
          WRITE(numout,*) ' '
          WRITE(numout,*) ' namelist for OPT model'
          WRITE(numout,*) ' ***********************'
          WRITE(numout,*) ' '
      ENDIF
C
      numnat=80
C
      clname='namelist.trc.sms'
      CALL ctlopn(numnat,clname,clold,clfor,clseq,
     $    ilseq,ilu,iused,numout,lwp,1)
C
C
C 1.5 namelist natopt : parameters for optic
C
      xkg0  = 0.
      xkr0  = 0.
      xkgp  = 0.
      xkrp  = 0.
      xlg   = 0.
      xlr   = 0.
      rpig  = 0.
C
      READ(numnat,natopt)
C
      IF(lwp) THEN
          WRITE(numout,*) 'natopt'
          WRITE(numout,*) ' '
          WRITE(numout,*) ' green   water absorption coeff  xkg0  = '
     $        ,xkg0
          WRITE(numout,*) ' red water absorption coeff      xkr0  = '
     $        ,xkr0
          WRITE(numout,*) ' pigment red absorption coeff    xkrp  = '
     $        ,xkrp
          WRITE(numout,*) ' pigment green absorption coeff  xkgp  = '
     $        ,xkgp
          WRITE(numout,*) ' green chl exposant              xlg   = '
     $        ,xlg
          WRITE(numout,*) ' red   chl exposant              xlr   = '
     $        ,xlr
          WRITE(numout,*) ' chla/chla+phea ratio            rpig  = '
     $        ,rpig
          WRITE(numout,*) ' '

      ENDIF
C
#if defined key_trc_diabio
C
C NAMELIST : natdbi 
C
C default name for biological trends : short and long name, units
C
      DO ji=1,jpdiabio
        IF (ji.lt.10) THEN 
            WRITE (ctrbio(ji),'("BIO_",I1)') ji
        ELSE IF (ji.lt.100) THEN
            WRITE (ctrbio(ji),'("BIO_",I2)') ji
        ELSE
            WRITE (ctrbio(ji),'("BIO_",I3)') ji
        ENDIF
        WRITE (ctrbil(ji),'("BIOLOGICAL TREND NUMBER ",I2)') ji
        ctrbiu(ji)='mmoleN/m3/s'
      END DO 
C
      nwritebio = 10
C
      READ(numnat,natdbi)
C
      IF(lwp) THEN
          WRITE(numout,*) 'natdbi'
          WRITE(numout,*) ' '
          WRITE(numout,*)
     $        ' frequency of outputs for biological outputs = '
     $        ,nwritebio
          WRITE(numout,*) ' '
          DO ji=1,jpdiabio
            WRITE(numout,*)
     $          'name of biological trend number :',ji,' : ',ctrbio(ji)  
            WRITE(numout,*) ctrbil(ji)  
            WRITE(numout,*) ' in unit = ',ctrbiu(ji)
          END DO 
      END IF 
C
#endif

#else
C
C no passive tracers
C
#endif
