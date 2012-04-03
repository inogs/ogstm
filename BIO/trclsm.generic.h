C $Id: trclsm.generic.h,v 1.2 2009-09-11 09:20:56 cvsogs01 Exp $
CCC
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

#if defined key_passivetrc && defined key_trc_npzd
      INTEGER ji
      INTEGER ilu,iused(1,100),ilseq
      CHARACTER*21 clold,clfor,clseq,clnew,cldir,clunf
      CHARACTER*32 clname
CC----------------------------------------------------------------------
CC statement functions
CC ===================

!                         #include "stafun.h"

CC
CCC---------------------------------------------------------------------
CCC  OPA8, LODYC (15/11/96)
CCC---------------------------------------------------------------------
C
C 0. initializations
C ------------------
C
      namelist/natbio/apmin,azmin,anmin,admin,
     $    redf,slopet,toptp,aknut,rcchl,
     $    rgamma,toptgz,tmaxgz,rgz,
     $    rppz,taus,aks,filmax,rpnaz,rdnaz,eggzoo,tauzn,
     $    tmmaxp,tmminp,tmmaxz,tmminz,anumin,afdmin,taudn,
     $    vsed,tmumax,aki,tmaxr,tminr,sedlam,sedlostpoc
      namelist/natopt/xkg0,xkr0,xkgp,xkrp,xlg,xlr,rpig
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
          WRITE(numout,*) ' namelist for generic model'
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
C
C 1.4 namelist natbio : biological parameters
C
C     apmin = 0
C     azmin = 0
C     anmin = 0
C     admin = 0
C     redf  = 0
C     slopet= 0
C     toptp = 0
C     aknut = 0
C     rcchl = 0
C     rgamma= 0
C     toptgz= 0
C     tmaxgz= 0
C     rgz   = 0
C     rppz  = 0
C     taus  = 0
C     aks   = 0
C     filmax= 0
C     rpnaz = 0
C     rdnaz = 0
C     eggzoo= 0
C     tauzn = 0
C     tmmaxp= 0
C     tmminp= 0
C     tmmaxz= 0
C     tmminz= 0
C     anumin= 0
C     afdmin= 0
C     taudn = 0
C     vsed  = 0
C     tmumax= 0
C     aki   = 0
C     tmaxr   = 1./(     4.*rday)*0.
C     tminr   = 1./(24.*30.*rday)*0.
C     sedlam=0
C     sedlostpoc=0

C
C
C     READ(numnat,natbio)
C
C     IF(lwp) THEN
C         WRITE(numout,*) 'natbio'
C         WRITE(numout,*) ' '
C         WRITE(numout,*)
C    $        ' minimum phytoplancton concentration  apmin =', apmin
C         WRITE(numout,*)
C    $        ' minimum zooplancton   concentration  azmin =', azmin
C         WRITE(numout,*)
C    $        ' minimum nutrients     concentration  anmin =', anmin
C         WRITE(numout,*)
C    $        ' minimum detritus      concentration  admin =', admin
C         WRITE(numout,*)
C    $        ' redfield ratio  c:n                   redf =', redf
C         WRITE(numout,*)
C    $        ' van t hoff coefficient              slopet =', slopet
C         WRITE(numout,*)
C    $        ' optimal photosynthesis temperature   toptp =', toptp
C         WRITE(numout,*)
C    $        ' half-saturation nutrient             aknut =', aknut
C         WRITE(numout,*)
C    $        ' carbone/chlorophyl ratio             rcchl =', rcchl
C         WRITE(numout,*)
C    $        ' phytoplankton exudation fraction    rgamma =', rgamma
C         WRITE(numout,*)
C    $        ' optimal temperature for zoo growth  toptgz =', toptgz
C         WRITE(numout,*)
C    $        ' maximal temperature for zoo growth  tmaxgz =', tmaxgz
C         WRITE(numout,*)
C    $        ' widtht of zoo temperature FUNCTION     rgz =', rgz
C         WRITE(numout,*)
C    $        ' zoo preference for phyto              rppz =', rppz
C         WRITE(numout,*)
C    $        ' maximal zoo grazing rate              taus =',86400*taus
C         WRITE(numout,*)
C    $        ' half saturation constant for zoo food  aks =', aks
C         WRITE(numout,*)
C    $        ' maximal mass clearance rate for zoo filmax =', filmax
C         WRITE(numout,*)
C    $        ' non-assimilated phyto by zoo         rpnaz =', rpnaz
C         WRITE(numout,*)
C    $        ' non-assimilated detritus by zoo      rdnaz =', rdnaz
C         WRITE(numout,*)
C    $        ' minimum  for zoo concentration      eggzoo =', eggzoo
C         WRITE(numout,*)
C    $        ' zoo specific excretion rate          tauzn =',86400
C    $        *tauzn
C         WRITE(numout,*)
C    $        ' maximal phyto mortality rate        tmmaxp =',86400
C    $        *tmmaxp
C         WRITE(numout,*)
C    $        ' minimal phyto mortality rate        tmminp =',86400
C    $        *tmminp
C         WRITE(numout,*)
C    $        ' maximal zoo mortality rate          tmmaxz =',86400
C    $        *tmmaxz
C         WRITE(numout,*)
C    $        ' minimal zoo mortality rate          tmminz =',86400
C    $        *tmminz
C         WRITE(numout,*)
C    $        ' nutrient threshold for phyto mort   anumin =', anumin
C         WRITE(numout,*)
C    $        ' food threshold for zoo mort         afdmin =', afdmin
C         WRITE(numout,*)
C    $        ' detrital breakdown rate              taudn =',86400
C    $        *taudn
C         WRITE(numout,*)
C    $        ' detritus sedimentation speed          vsed =',86400*vsed
C         WRITE(numout,*)
C    $        ' phyto max growth rate               tmumax =',86400
C    $        *tmumax
C         WRITE(numout,*)
C    $        ' light hlaf saturation constant         aki =', aki
C         WRITE(numout,*)
C    $        ' maximum damping for d z or p         tmaxr =', tmaxr
C         WRITE(numout,*)
C    $        ' minimum damping for d z or p         tminr =', tminr
C         WRITE(numout,*)
C    $        ' time coeff of POC in sediments      sedlam =', sedlam
C         WRITE(numout,*)
C    $        ' Sediment geol loss for POC  sedlostpoc =', sedlostpoc
C     ENDIF
C
C 1.5 namelist natopt : parameters for optic
C
C     xkg0  = 0.
C     xkr0  = 0.
C     xkgp  = 0.
C     xkrp  = 0.
C     xlg   = 0.
C     xlr   = 0.
C     rpig  = 0.
C
C     READ(numnat,natopt)
C
C     IF(lwp) THEN
C         WRITE(numout,*) 'natopt'
C         WRITE(numout,*) ' '
C         WRITE(numout,*) ' green   water absorption coeff  xkg0  = '
C    $        ,xkg0
C         WRITE(numout,*) ' red water absorption coeff      xkr0  = '
C    $        ,xkr0
C         WRITE(numout,*) ' pigment red absorption coeff    xkrp  = '
C    $        ,xkrp
C         WRITE(numout,*) ' pigment green absorption coeff  xkgp  = '
C    $        ,xkgp
C         WRITE(numout,*) ' green chl exposant              xlg   = '
C    $        ,xlg
C         WRITE(numout,*) ' red   chl exposant              xlr   = '
C    $        ,xlr
C         WRITE(numout,*) ' chla/chla+phea ratio            rpig  = '
C    $        ,rpig
C         WRITE(numout,*) ' '

C     ENDIF
C
#if defined key_trc_diabio
C
C NAMELIST : natdbi 
C
C default name for biological trends : short and long name, units
C
C     DO ji=1,jpdiabio
C       IF (ji.lt.10) THEN 
C           WRITE (ctrbio(ji),'("BIO_",I1)') ji
C       ELSE IF (ji.lt.100) THEN
C           WRITE (ctrbio(ji),'("BIO_",I2)') ji
C       ELSE
C           WRITE (ctrbio(ji),'("BIO_",I3)') ji
C       ENDIF
C       WRITE (ctrbil(ji),'("BIOLOGICAL TREND NUMBER ",I2)') ji
C       ctrbiu(ji)='mmoleN/m3/s'
C     END DO 
C
C     nwritebio = 10
C
C     READ(numnat,natdbi)
C
C     IF(lwp) THEN
C         WRITE(numout,*) 'natdbi'
C         WRITE(numout,*) ' '
C         WRITE(numout,*)
C    $        ' frequency of outputs for biological outputs = '
C    $        ,nwritebio
C         WRITE(numout,*) ' '
C         DO ji=1,jpdiabio
C           WRITE(numout,*)
C    $          'name of biological trend number :',ji,' : ',ctrbio(ji)  
C           WRITE(numout,*) ctrbil(ji)  
C           WRITE(numout,*) ' in unit = ',ctrbiu(ji)
C         END DO 
C     END IF 
C
Cendif

#else
C
C no passive tracers
C
#endif
