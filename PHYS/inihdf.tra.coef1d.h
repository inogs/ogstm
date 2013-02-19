CC $Header: /cvsroot/opatm-bfm/opa_model/OPA/inihdf.tra.coef1d.h,v 1.2 2009-09-11 09:20:56 cvsogs01 Exp $
CCC
CCC             inihdf.tra.coef1d.h
CCC           ***********************
CCC
CCC   defined key :  'key_trahdfcoef1d'
CCC   ============
CCC
CC      1D eddy diffusivity coefficients ( depth )
CC      --------------------------------
CC
CC       biharmonic operator    : ahtt = defined at T-level
CC                                ahtu,ahtv,ahtw never used
CC
CC       harmonic operator      : ahtt never used
CC          -1-  iso-model level: ahtu = ahtv defined at T-level
CC                                ahtw never used
CC          -2-  isopycnal or   : ahtu = ahtv defined at T-level
CC               geopotential     ahtw defined at w-level
CC
CC       eddy induced velocity
CC         always harmonic      : aeiu = aeiv defined at T-level
CC				  aeiw defined at w-level
CC
CCC---------------------------------------------------------------------
CCC  OPA8, LODYC (1997)
CCC---------------------------------------------------------------------
C
      IF(lwp)WRITE(numout,*)
#ifdef key_trahdfeiv
      IF(lwp)WRITE(numout,*) ' inihdf: 1D eddy diffusivity and eddy',
     $                       ' induced velocity coefficients'
#  else
      IF(lwp)WRITE(numout,*) ' inihdf: 1D eddy diffusivity coefficient '
#endif
      IF(lwp)WRITE(numout,*) ' ======  --'
      IF(lwp)WRITE(numout,*)
C    
C ... initialization of the profile
C
C   ... ahts, ahtf: surface and bottom values
      zahts = 2000.
      zahtf =  500.
CCC 22/3/2004 Paolo dynamical models Mediterranean Scale Settings
       zahts = aht0
       zahtf = aht0/2.
CCC 
C   ... zkah, zahr: depth of the inflection pt and width of inflection
CCC 22/3/2004 Paolo dynamical models Mediterranean Scale Settings
CCC       zkah =  -300.
CCC       zahr =   300.
          zkah =  -150.
          zahr =   150.
CCC
C   ... computation coefficients
      za00 = tanh( (-fsdept(1,1,1  )-zkah) / zahr )
      za01 = tanh( (-fsdept(1,1,jpk)-zkah) / zahr )
      zahf = ( zahts-zahtf ) / ( za00 - za01 )
      zahs = zahts - zahf * za00
C
#if defined key_trahdfbilap
C
C biharmonic operator : 
C ==================== 
C
C ... set ahtt at T-level
      DO jk = 1, jpk
CC-CC        ahtt(jk) = zahs + zahf * tanh( (-fsdept(1,1,jk)-zkah) / zahr )
        ahtt(jk) = aht0
      END DO
C
C ... control print
      IF(lwp) THEN
          WRITE(numout,*)
          WRITE(numout,*) '         aht profile at T-level : '
          WRITE(numout,*)
          WRITE(numout,9200)
          DO jk = 1, jpk
            WRITE(numout,9210) jk, ahtt(jk), fsdept(1,1,jk)
          END DO
      ENDIF
 9200 FORMAT(' level    aht          depth t-level ' )
 9210 FORMAT(i6,2f12.4)
C
#  else
C
C harmonic operator :
C ================== 
C
C ... set ahtu = ahtv at T-level, and ahtw at w-level
      DO jk = 1, jpk
        ahtu(jk) = zahs + zahf * tanh( (-fsdept(1,1,jk)-zkah) / zahr )
        ahtv(jk) = ahtu(jk)
        ahtw(jk) = zahs + zahf * tanh( (-fsdepw(1,1,jk)-zkah) / zahr )
      END DO
C
C ... control print
      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,*) '         aht profile at T-level : '
      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,9200)
      DO jk = 1, jpk
        IF(lwp)WRITE(numout,9210) jk, ahtu(jk), fsdept(1,1,jk)
      END DO
 9200 FORMAT(' level    aht          depth t-level ' )
 9210 FORMAT(i6,2f12.4)
      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,*) '         aht profile at W-level : '
      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,9220)
      DO jk = 1, jpk
        IF(lwp)WRITE(numout,9210) jk, ahtw(jk), fsdepw(1,1,jk)
      END DO
 9220 FORMAT('  jk      aht          depth w-level ' )
C
#    ifdef key_trahdfeiv
C ... set aeiu = aeiv and set aeiw (here same profile as on aht)
C
      DO jk = 1, jpk
        aeiu(jk) = ahtu(jk)
        aeiv(jk) = aeiu(jk)
        aeiw(jk) = ahtw(jk)
      END DO
C
C ... Control print
      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,*) '         aeiv profile at T-level : '
      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,9230)
      DO jk = 1, jpk
        IF(lwp)WRITE(numout,9240) jk, aeiu(jk), fsdept(1,1,jk)
      END DO
 9230 FORMAT('  jk      aeiv           depth t-level ' )
 9240 FORMAT(i6,2f12.4)
#    endif
C
#endif