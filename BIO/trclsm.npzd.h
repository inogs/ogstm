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
      CHARACTER(LEN=21) clold,clfor,clseq,clnew,cldir,clunf
      CHARACTER(LEN=32) clname

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
          WRITE(numout,*) ' namelist for NPZD model'
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
      apmin = 0
      azmin = 0
      anmin = 0
      admin = 0
      redf  = 0
      slopet= 0
      toptp = 0
      aknut = 0
      rcchl = 0
      rgamma= 0
      toptgz= 0
      tmaxgz= 0
      rgz   = 0
      rppz  = 0
      taus  = 0
      aks   = 0
      filmax= 0
      rpnaz = 0
      rdnaz = 0
      eggzoo= 0
      tauzn = 0
      tmmaxp= 0
      tmminp= 0
      tmmaxz= 0
      tmminz= 0
      anumin= 0
      afdmin= 0
      taudn = 0
      vsed  = 0
      tmumax= 0
      aki   = 0
      tmaxr   = 1./(     4.*rday)*0.
      tminr   = 1./(24.*30.*rday)*0.
      sedlam=0
      sedlostpoc=0

C
C
      READ(numnat,natbio)
C
      IF(lwp) THEN
          WRITE(numout,*) 'natbio'
          WRITE(numout,*) ' '
          WRITE(numout,*)
     $        ' minimum phytoplancton concentration  apmin =', apmin
          WRITE(numout,*)
     $        ' minimum zooplancton   concentration  azmin =', azmin
          WRITE(numout,*)
     $        ' minimum nutrients     concentration  anmin =', anmin
          WRITE(numout,*)
     $        ' minimum detritus      concentration  admin =', admin
          WRITE(numout,*)
     $        ' redfield ratio  c:n                   redf =', redf
          WRITE(numout,*)
     $        ' van t hoff coefficient              slopet =', slopet
          WRITE(numout,*)
     $        ' optimal photosynthesis temperature   toptp =', toptp
          WRITE(numout,*)
     $        ' half-saturation nutrient             aknut =', aknut
          WRITE(numout,*)
     $        ' carbone/chlorophyl ratio             rcchl =', rcchl
          WRITE(numout,*)
     $        ' phytoplankton exudation fraction    rgamma =', rgamma
          WRITE(numout,*)
     $        ' optimal temperature for zoo growth  toptgz =', toptgz
          WRITE(numout,*)
     $        ' maximal temperature for zoo growth  tmaxgz =', tmaxgz
          WRITE(numout,*)
     $        ' widtht of zoo temperature FUNCTION     rgz =', rgz
          WRITE(numout,*)
     $        ' zoo preference for phyto              rppz =', rppz
          WRITE(numout,*)
     $        ' maximal zoo grazing rate              taus =',86400*taus
          WRITE(numout,*)
     $        ' half saturation constant for zoo food  aks =', aks
          WRITE(numout,*)
     $        ' maximal mass clearance rate for zoo filmax =', filmax
          WRITE(numout,*)
     $        ' non-assimilated phyto by zoo         rpnaz =', rpnaz
          WRITE(numout,*)
     $        ' non-assimilated detritus by zoo      rdnaz =', rdnaz
          WRITE(numout,*)
     $        ' minimum  for zoo concentration      eggzoo =', eggzoo
          WRITE(numout,*)
     $        ' zoo specific excretion rate          tauzn =',86400
     $        *tauzn
          WRITE(numout,*)
     $        ' maximal phyto mortality rate        tmmaxp =',86400
     $        *tmmaxp
          WRITE(numout,*)
     $        ' minimal phyto mortality rate        tmminp =',86400
     $        *tmminp
          WRITE(numout,*)
     $        ' maximal zoo mortality rate          tmmaxz =',86400
     $        *tmmaxz
          WRITE(numout,*)
     $        ' minimal zoo mortality rate          tmminz =',86400
     $        *tmminz
          WRITE(numout,*)
     $        ' nutrient threshold for phyto mort   anumin =', anumin
          WRITE(numout,*)
     $        ' food threshold for zoo mort         afdmin =', afdmin
          WRITE(numout,*)
     $        ' detrital breakdown rate              taudn =',86400
     $        *taudn
          WRITE(numout,*)
     $        ' detritus sedimentation speed          vsed =',86400*vsed
          WRITE(numout,*)
     $        ' phyto max growth rate               tmumax =',86400
     $        *tmumax
          WRITE(numout,*)
     $        ' light hlaf saturation constant         aki =', aki
          WRITE(numout,*)
     $        ' maximum damping for d z or p         tmaxr =', tmaxr
          WRITE(numout,*)
     $        ' minimum damping for d z or p         tminr =', tminr
          WRITE(numout,*)
     $        ' time coeff of POC in sediments      sedlam =', sedlam
          WRITE(numout,*)
     $        ' Sediment geol loss for POC  sedlostpoc =', sedlostpoc
      ENDIF
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
