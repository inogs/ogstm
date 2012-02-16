C $Id: p3zbio.F,v 1.2 2009-09-11 09:20:56 cvsogs01 Exp $
      SUBROUTINE p3zbio
#if defined key_passivetrc && defined key_trc_p3zd
CCC   ------------------------------------------------------------------
CCC   
CCC   ROUTINE p3zbio
CCC   ******************
CCC   
CCC   
CC
CC     PURPOSE.
CC     --------
CC          *P3ZBIO* MODELS PRODUCTION OF BIOGENIC MATTER (POC '' SOFT
CC                   TISSUE'' AND CACO3 PARTICLES ''HARD PARTS'')
CC                   AND ITS DISTRIBUTION IN WATER COLUMN
CC
CC     METHOD.
CC     -------
CC          IN THE SURFACE LAYER POC IS PRODUCED ACCORDING TO
CC     NURTRIENTS AVAILABLE AND GROWTH CONDITIONS. NUTRIENT UPTAKE
CC     KINETICS FOLLOW MICHAELIS-MENTON FORMULATION. PROPORTIONAL
CC     TO THE AMOUNT OF ORGANIC MATTER, CACO3 HARD PARTS ARE PRODUCED.
CC     THE TOTAL PARTICLE AMOUNT PRODUCED, IS DISTRIBUTED IN THE WATER
CC     COLUMN BELOW THE SURFACE LAYER.
CC
CC     EXTERNALS.
CC     ----------
CC          NONE.
CC
CC     REFERENCE.
CC     ----------
CC
CC          BACASTOW, R., AND E. MAIER-REIMER (1985)
CC          CIRCULATION MODEL OF THE OCEAN CARBON CYCLE.
CC          1. DESCRIPTION OF THE MODEL, PP. 224-232.
CC          2. COMPARISON OF THE MODEL RESULTS WITH OBSERVATIONAL DATA,
CC          PP. 233-240.
CC          IN: "ATMOSPHERIC CARBON DIOXIDE - ITS SOURCES, SINKS, AND
CC          GLOBAL tranSPORT", KANDERSTEG, 2 TO 6 SEPTEMBER 1985,
CC          COMMISSION ON ATMOSPHERIC CHEMISTRY AND GLOBAL POLLUTION,
CC          INTERNATIONAL ASSOCIATION OF METEOROLOGY AND ATMOSPHERIC PHYSICS.
CC
CC          DUGDALE. R.C. (1967)
CC          NUTRIENT LIMITATION IN THE SEA: DYNAMICS, IDENTIFICATION
CC          AND SIGNIFICANCE.
CC          LIMNOLOGY AND OCEANOGRAPHY, VOL.12, 685-695.
CC
CC          PARSONS, T.R., AND M. TAKAHASHI (1973)
CC          BIOLOGICAL OCEANOGRAPHIC PROCESSES.
CC          PERGAMON PRESS, 186 PP.
CC
CC   MODIFICATIONS:
CC   --------------
CC      original      : 1998    O. Aumont
CC      modifications : 1999    C. Le Quere
CC      modifications : 1999    O. Aumont
CC      modifications : 2001    O. Aumont
CC ----------------------------------------------------------------
CC parameters and commons
CC ======================
C     DIR$ NOLIST
c#include "parameter.h"
c#include "common.h"
       USE myalloc
C     DIR$ LIST
CC+CC  Implicit typing is never allowed
        IMPLICIT NONE
CC+CC  Implicit typing is never allowed

CC-----------------------------------------------------------------
c -----
CC local declarations
CC ==================
C     
      INTEGER ji, jj, jk
      INTEGER kmin(jpi,jpj)

      REAL(8) silpot, expofa
      REAL(8)  calpot, silfra
      REAL(8) remip,remik
C     
      REAL(8) emoy(jpi,jpj),dipnmoy(jpi,jpj),zmeu(jpi,jpj)
      REAL(8) orem(jpi,jpj,jpk),olimi(jpi,jpj,jpk),intpz(jpi,jpj)
      REAL(8) phosph2,zoo2,oxygen2
      REAL(8) phyto2,poc212
      REAL(8) compaph(jpi,jpj,jpk),compaz
C     
      REAL(8) parlux
      REAL(8) graze,prefc,prefp
C     
C     SET HALF PRECISION CONSTANTS
C-----------------------------
C     
      e1       = 0.
      e2       = 0.
      zero     = 0.
      one      = 1.
      two      = 2.
      parlux   = 0.225
      
      intpz       = 0.   
C     
C     Initialisation of variables used to compute PAR
C     -----------------------------------------------
C     
      emoy        = 0.
      dipnmoy     = 0.
      etot        = 0.
C     
C     Light penetration in the water column
C     -------------------------------------
      DO jj = 1,jpj
        DO ji = 1,jpi
          e1(ji,jj,1)=parlux*qsr(ji,jj)
          e2(ji,jj,1)=parlux*qsr(ji,jj)
        END DO
      END DO

      DO jk = 1,16
        DO jj = 1,jpj
          DO ji = 1,jpi
            IF (tmask(ji,jj,jk).NE.0) THEN
C     
C     Separation in two light bands: red and green
C     --------------------------------------------
C     
                e1(ji,jj,jk+1) = e1(ji,jj,jk)*
     &              exp(-(ekw1+ekc*trn(ji,jj,jk,jpphy)*12e6
     &              /dipn(ji,jj,jk))*e3t(jk)/2.)
                e2(ji,jj,jk+1) = e2(ji,jj,jk)*
     &              exp(-(ekw2+ekc*trn(ji,jj,jk,jpphy)*12e6
     &              /dipn(ji,jj,jk))*e3t(jk)/2.)
                etot(ji,jj,jk) = e1(ji,jj,jk+1)+e2(ji,jj,jk+1)
C     
C     Computation of irradiance below level T
C     ---------------------------------------
C     
                e1(ji,jj,jk+1)=e1(ji,jj,jk+1)*
     $              exp(-(ekw1+ekc*trn(ji,jj,jk,jpphy)*12e6
     $              /dipn(ji,jj,jk))*e3t(jk)/2.)
                e2(ji,jj,jk+1)=e2(ji,jj,jk+1)*
     $              exp(-(ekw2+ekc*trn(ji,jj,jk,jpphy)*12e6
     $              /dipn(ji,jj,jk))*e3t(jk)/2.)
C     
C     Computation of C/Chl ratio (doney et al., 1996)
C     -----------------------------------------------
C     
                dipn(ji,jj,jk)=trn(ji,jj,jk,jppo4)
     $              /(trn(ji,jj,jk,jppo4)+conc0)
     $              *(2.5-1.5*MIN(etot(ji,jj,jk)/90.,1.))/7.6/12.
     $              +1.E-30
                dipn(ji,jj,jk)=MIN(150.,1./dipn(ji,jj,jk))  
            ENDIF
          END DO
        END DO
      END do 
C     
C     Initialisation of the euphotic depth
C     
      DO jj = 1,jpj
        DO ji = 1,jpi
          kmin(ji,jj)=1
          zmeu(ji,jj)=gdept(16)
        END DO
      END DO     
C     
C     Computation of the euphotic depth
C     ---------------------------------
C     
      DO jk = 2,16
        DO jj = 1,jpj
          DO ji = 1,jpi
            IF (etot(ji,jj,jk).GE.0.0045*qsr(ji,jj)) THEN 
                zmeu(ji,jj) = gdept(jk)
                IF (gdept(jk).LE.hmld(ji,jj)) THEN
                    kmin(ji,jj)=jk
                ENDIF
            ENDIF
          END DO
        END DO
      END DO


      DO jk = 1,16
        DO jj = 1,jpj
          DO ji = 1,jpi
            IF (jk.LE.kmin(ji,jj)) THEN
C     
C     C/CHL ET PAR MOYENS SUR LA COUCHE DE MELANGE
C     
                emoy(ji,jj)=emoy(ji,jj)+etot(ji,jj,jk)/
     $              (float(kmin(ji,jj))+1.E-15)
                dipnmoy(ji,jj)=dipnmoy(ji,jj)+dipn(ji,jj,jk)/
     $              (float(kmin(ji,jj))+1.E-15)
            ENDIF
          END DO
        END DO
      END DO


      DO jk = 1,jpkb
        DO jj = 1,jpj
          DO ji = 1,jpi
            IF (jk.LE.kmin(ji,jj)) THEN
                etot(ji,jj,jk)=emoy(ji,jj)
                dipn(ji,jj,jk)=dipnmoy(ji,jj)
            ENDIF
C     
#    if defined key_trc_dia3d
            trc3d(ji,jj,jk,1) = dipn(ji,jj,jk)
#    endif
          END DO
        END DO
      END DO

#    if defined key_diatrdtrc
      DO jj = 1,jpj
        DO ji = 1,jpi
          trc2d(ji,jj,11) = dipnmoy(ji,jj)
          trc2d(ji,jj,12) = emoy(ji,jj)
        END DO
      END DO
#    endif

C     
      DO  jk = 1,jpkb
        DO  jj = 1,jpj
          DO  ji = 1,jpi
            IF (tmask(ji,jj,jk).NE.0) THEN
C     
C     Computation of phyto development constant
C     -----------------------------------------
C     
                prbio(ji,jj,jk) = 0.6/rjjss*(1.066)**(tn(ji,jj,jk))
C     
C     Computation of production function
C     ----------------------------------
C     
         prbio(ji,jj,jk) = prbio(ji,jj,jk)*(1.-exp(-pislope/rjjss*
     &     etot(ji,jj,jk)/prbio(ji,jj,jk)))*exp(-betslope/rjjss*
     &     etot(ji,jj,jk)/prbio(ji,jj,jk))
#    if defined key_off_degrad
     &              *facvol(ji,jj,jk)
#    endif
C     
C     Mixed-layer effect on production
C     --------------------------------
C     
                IF (hmld(ji,jj).ge.2*zmeu(ji,jj)) THEN
                    prbio(ji,jj,jk) = prbio(ji,jj,jk)*0.5
                    
                ELSEIF (hmld(ji,jj).le.zmeu(ji,jj)) THEN
                    prbio(ji,jj,jk) = prbio(ji,jj,jk)
                ELSE
                    prbio(ji,jj,jk) = prbio(ji,jj,jk)*(1.-
     &                  0.5*(hmld(ji,jj)/zmeu(ji,jj)-1))
                ENDIF
            ENDIF
          ENDDO
        ENDDO

C     
        DO jj = 1,jpj
          DO ji = 1,jpi
            IF (tmask(ji,jj,jk).ne.0) THEN
C     
C     Exsudation of Zoo towards DOC
C     -----------------------------
C     
                compaz = max((trn(ji,jj,jk,jpzoo)-0.01E-6),0.)
                compaph(ji,jj,jk) = max((trn(ji,jj,jk,jpphy)-0.01E-6),
     $              0.)
                respz(ji,jj,jk) = resrat/rjjss*trn(ji,jj,jk,jpzoo)/
     &              (1.E-6+trn(ji,jj,jk,jpzoo))*rfact*compaz*
     $              tmask(ji,jj,jk)
#    if defined key_off_degrad
     &              *facvol(ji,jj,jk)
#    endif
c     
C     
C     Computation of the fast remineralised fraction as a function of nutrients
C     -------------------------------------------------------------------------
C     
                eps1(ji,jj,jk) =epsbio+0.*7.6e-6/(7.6e-6+
     $              trn(ji,jj,jk,jppo4))
C     
C     Squared mortality of Phyto similar to a sedimentation term during
c blooms
C     (Doney et al. 1996)
C     
                respp(ji,jj,jk) = wchl/rjjss*1e6*
     $              trn(ji,jj,jk,jpphy)**2*rfact
     &              *tmask(ji,jj,jk)
#    if defined key_off_degrad
     &              *facvol(ji,jj,jk)
#    endif
C     
C     Phytoplankton mortality
C     -----------------------
C     
                tortp(ji,jj,jk) = mprat/rjjss*trn(ji,jj,jk,jpphy)
     $              /(1.e-6+trn(ji,jj,jk,jpphy))
     &              *compaph(ji,jj,jk)*rfact*tmask(ji,jj,jk)
#    if defined key_off_degrad
     &              *facvol(ji,jj,jk)
#    endif
C     
C     Zooplankton mortality
C     ---------------------
C     
                tortz(ji,jj,jk) = mzrat/rjjss*trn(ji,jj,jk,jpzoo)
     $              /(1.E-6+trn(ji,jj,jk,jpzoo))
     &              *compaz*rfact*tmask(ji,jj,jk)
#    if defined key_off_degrad
     &              *facvol(ji,jj,jk)
#    endif
C     
C     Zooplankton grazing
C     -------------------
C     
                graze = grazrat/rjjss*rfact*tmask(ji,jj,jk)
#    if defined key_off_degrad
     &              *facvol(ji,jj,jk)
#    endif
C     
C     Preference of zooplankton for Phyto and POC (Fasham et al. 1990)
C     
                prefc = xprefc*trn(ji,jj,jk,jpphy)/
     $              (xprefc*trn(ji,jj,jk,jpphy)
     &              +xprefp*trn(ji,jj,jk,jppoc)+1E-15)

                prefp = xprefp*trn(ji,jj,jk,jppoc)/
     $              (xprefc*trn(ji,jj,jk,jpphy)
     &              +xprefp*trn(ji,jj,jk,jppoc)+1E-15)

                grazp(ji,jj,jk) = graze*prefc*compaph(ji,jj,jk)/
     $              (xkgraz+prefc*trn(ji,jj,jk,jpphy)+prefp*
     &              trn(ji,jj,jk,jppoc))*trn(ji,jj,jk,jpzoo)
                

                grazpoc(ji,jj,jk) = graze*prefp*trn(ji,jj,jk,jppoc)/
     &              (xkgraz+prefc*trn(ji,jj,jk,jpphy)+prefp*
     &              trn(ji,jj,jk,jppoc))*trn(ji,jj,jk,jpzoo)
C     
C     
C     Sedimentation of Phyto and POC
C     ------------------------------
C     
                sinking(ji,jj,jk+1) = wsbio/rjjss*trn(ji,jj,jk,jppoc)
     &              *rfact*tmask(ji,jj,jk+1)
#    if defined key_off_degrad
     &              *facvol(ji,jj,jk)
#    endif
                nu(ji,jj,jk+1) = smax/rjjss*(1.-tanh(conc0*0.0481e12*
     &              trn(ji,jj,jk,jppo4)))*trn(ji,jj,jk,jpphy)*rfact*
     &              tmask(ji,jj,jk+1)
#    if defined key_off_degrad
     &              *facvol(ji,jj,jk)
#    endif
C     
C     Remineralization of DOC and POC
C     
                remik=1.64*SPOCRI*(1.-tmask(ji,jj,jk+1))+xremik
                remik=remik/rjjss*rfact*tmask(ji,jj,jk)
     &              *trn(ji,jj,jk,jppo4)/(trn(ji,jj,jk,jppo4)+32.E-6)
     &              *trn(ji,jj,jk,jpdoc)/(trn(ji,jj,jk,jpdoc)+15.E-6)
#    if defined key_off_degrad
     &              *facvol(ji,jj,jk)
#    endif
                remip=xremip/rjjss*rfact*tmask(ji,jj,jk)
#    if defined key_off_degrad
     &              *facvol(ji,jj,jk)
#    endif
                olimi(ji,jj,jk)=remik*trn(ji,jj,jk,jpdoc) 
                orem(ji,jj,jk)=remip*trn(ji,jj,jk,jppoc)
            ENDIF
          END DO
        END DO
      END DO
C     
      DO jk = 1,jpkb
        DO jj = 1,jpj
          DO ji = 1,jpi
            IF (tmask(ji,jj,jk).NE.0) THEN
C     
C     Computation of Primary Production 
C     ---------------------------------
C     
           prorca(ji,jj,jk) = prbio(ji,jj,jk)*trn(ji,jj,jk,jppo4)*
     &              compaph(ji,jj,jk)/(conc0+trn(ji,jj,jk,jppo4))*
     &              rfact
C     
C     Evolution of PO4
C     ----------------
C     
                phosph2 = trn(ji,jj,jk,jppo4)-prorca(ji,jj,jk)
     &              +eps1(ji,jj,jk)*tortz(ji,jj,jk)
     &              +olimi(ji,jj,jk)
C     
C     Nullity test for PO4
C     --------------------
C     
                prorca(ji,jj,jk)=prorca(ji,jj,jk)
     &              *(0.5+sign(0.5,phosph2))

C     
C     Evolution of Phyto
C     ------------------
C     
                phyto2 = trn(ji,jj,jk,jpphy)+prorca(ji,jj,jk)*
     &              (1.-excret)-(tortp(ji,jj,jk)+
     &              respp(ji,jj,jk))-grazp(ji,jj,jk)+
     &              (nu(ji,jj,jk)-nu(ji,jj,jk+1))/e3t(jk)
C     
C     Nullity test for Phyto
C     ----------------------
C     
                tortp(ji,jj,jk)=tortp(ji,jj,jk)
     $              *(0.5+sign(0.5,phyto2))
                respp(ji,jj,jk)=respp(ji,jj,jk)
     $              *(0.5+sign(0.5,phyto2))
                grazp(ji,jj,jk)=grazp(ji,jj,jk)
     $              *(0.5+sign(0.5,phyto2))
                prorca(ji,jj,jk)=prorca(ji,jj,jk)
     $              *(0.5+sign(0.5,phyto2))
                nu(ji,jj,jk+1)=nu(ji,jj,jk+1)
     $              *(0.5+sign(0.5,phyto2))
C     
C     Evolution of detritus
C     ---------------------
C     
                poc212 = trn(ji,jj,jk,jppoc)-grazpoc(ji,jj,jk)+unass*
     &              (grazp(ji,jj,jk)+grazpoc(ji,jj,jk))+
     &              (1.-eps1(ji,jj,jk))*tortz(ji,jj,jk)+
     &              tortp(ji,jj,jk)+respp(ji,jj,jk) +
     &              (sinking(ji,jj,jk)-sinking(ji,jj,jk+1))/e3t(jk)
     &              -orem(ji,jj,jk)
C     
C     Nullity test for POC
C     --------------------
C     
                grazpoc(ji,jj,jk)=grazpoc(ji,jj,jk)
     &              *(0.5+sign(0.5,poc212))
                sinking(ji,jj,jk+1)=sinking(ji,jj,jk+1)
     &              *(0.5+sign(0.5,poc212))
                orem(ji,jj,jk)=orem(ji,jj,jk) 
     &              *(0.5+sign(0.5,poc212))
C     
C     Evolution of Zooplankton
C     ------------------------
C     
                zoo2 = trn(ji,jj,jk,jpzoo)+
     &              (1-unass)*(grazp(ji,jj,jk)+grazpoc(ji,jj,jk))
     &              -(tortz(ji,jj,jk)+respz(ji,jj,jk))
C     
C     
C     Nullity test for Zooplankton
C     ----------------------------
C     
                tortz(ji,jj,jk)=tortz(ji,jj,jk)
     &              *(0.5+sign(0.5,zoo2))
                respz(ji,jj,jk)=respz(ji,jj,jk)
     &              *(0.5+sign(0.5,zoo2))
                grazp(ji,jj,jk)=grazp(ji,jj,jk)
     &              *(0.5+sign(0.5,zoo2))
                grazpoc(ji,jj,jk)=grazpoc(ji,jj,jk)
     &              *(0.5+sign(0.5,zoo2))
C     
C     Evolution of O2
C     ---------------
C     
                oxygen2 = trn(ji,jj,jk,jpoxy)+o2ut*
     &              (prorca(ji,jj,jk)-eps1(ji,jj,jk)*tortz(ji,jj,jk)
     &              -olimi(ji,jj,jk))
C     
                tortz(ji,jj,jk)=tortz(ji,jj,jk)
     &              *(0.5+sign(0.5,oxygen2))
                olimi(ji,jj,jk)=olimi(ji,jj,jk)
     &              *(0.5+sign(0.5,oxygen2))
C     
            ENDIF
          END DO
        END DO
      END DO
C     
C     Determination of tracers concentration as a function of 
C     biological sources and sinks
C     --------------------------------------------------------
C     
      DO jk = 1,jpkb
        DO jj = 1,jpj
          DO ji = 1,jpi
C     
            IF (tmask(ji,jj,jk).NE.0) THEN
C     
C     Evolution of PO4
C     ----------------
C     
                trn(ji,jj,jk,jppo4) = trn(ji,jj,jk,jppo4)-
     &              prorca(ji,jj,jk)+olimi(ji,jj,jk)+
     &              eps1(ji,jj,jk)*tortz(ji,jj,jk)
                prodt(ji,jj,jk) = prorca(ji,jj,jk)
     &              -eps1(ji,jj,jk)*tortz(ji,jj,jk)
     &              -olimi(ji,jj,jk)
C     
C     Evolution of Phytoplankton
C     --------------------------
C     
C     
                trn(ji,jj,jk,jpphy) = trn(ji,jj,jk,jpphy)+
     &              prorca(ji,jj,jk)*(1.-excret)-
     &              (tortp(ji,jj,jk)+respp(ji,jj,jk))-
     &              grazp(ji,jj,jk)+
     &              (nu(ji,jj,jk)-nu(ji,jj,jk+1))/e3t(jk)
C     
C     Evolution of Zooplankton
C     ------------------------
C     
                trn(ji,jj,jk,jpzoo) = trn(ji,jj,jk,jpzoo)+
     &              (1-unass)*(grazp(ji,jj,jk)+grazpoc(ji,jj,jk))
     &              -(tortz(ji,jj,jk)+respz(ji,jj,jk))
C     
C     Evolution of DOC
C     ----------------
C     
                trn(ji,jj,jk,jpdoc) = trn(ji,jj,jk,jpdoc)
     &              +respz(ji,jj,jk)+orem(ji,jj,jk)
     &              +excret*prorca(ji,jj,jk)-olimi(ji,jj,jk)
C     
C     Evolution of Detritus
C     ---------------------
C     
C     
                trn(ji,jj,jk,jppoc) = trn(ji,jj,jk,jppoc)-
     &              grazpoc(ji,jj,jk)+unass*
     &              (grazp(ji,jj,jk)+grazpoc(ji,jj,jk))+
     &              (1.-eps1(ji,jj,jk))*tortz(ji,jj,jk)+
     &              tortp(ji,jj,jk)+respp(ji,jj,jk)
     &              +(sinking(ji,jj,jk)-sinking(ji,jj,jk+1))/e3t(jk)
     &              -orem(ji,jj,jk)
C     
C     Evolution of O2
C     ---------------
C     
                trn(ji,jj,jk,jpoxy)=trn(ji,jj,jk,jpoxy)+
     &              o2ut*(prorca(ji,jj,jk)-olimi(ji,jj,jk)-
     &              eps1(ji,jj,jk)*tortz(ji,jj,jk))
C     
C     Vertical integral of phyto and zooplanton concentrations
C     used to compute calcite production
C     
                intpz(ji,jj)=intpz(ji,jj)+(trn(ji,jj,jk,jpphy)
     &              +trn(ji,jj,jk,jpzoo))*e3t(jk)*tmask(ji,jj,jk)
            ENDIF
          END DO
        END DO
      END DO
      
C     
C     Evolution of calcite and silicates as a function of the two tracers
C     -------------------------------------------------------------------
C     
      DO  jk = 1,jpkb
        DO  jj = 1,jpj
          DO  ji = 1,jpi
            IF (tmask(ji,jj,jk).ne.0) THEN
C     
C     potential production of calcite and biogenic silicate
C     ------------------------------------------------------
C     
                silpot = grosip*trn(ji,jj,jk,jpsil)/
     &              (trn(ji,jj,jk,jpsil)+xksi1)
                expofa=EXP(0.1*tn(ji,jj,jk)-10.)
                calpot=expofa/(1.+expofa)
                calfra=caco3r*calpot
                silfra=silpot*(1.-calpot)/0.75
C     
C     in situ production of calcite and biogenic silicate
C     ----------------------------------------------------
C     
                prcaca(ji,jj,jk)=calfra*(sinking(ji,jj,11)
     &              +nu(ji,jj,11))*(trn(ji,jj,jk,jpphy)+
     &              trn(ji,jj,jk,jpzoo))*tmask(ji,jj,1)
     &              /(intpz(ji,jj)+1.E-15)
                silpro(ji,jj,jk)=sicmax*silfra*
     &              (unass*(grazp(ji,jj,jk)+grazpoc(ji,jj,jk))
     &              +tortp(ji,jj,jk)+respp(ji,jj,jk))
C     
C     Compute variable Si/C ratio. Very simple formulation 
C     based on the assumption that this ratio decreased with
C     decreasing in situ Si concentrations
C     --------------------------------------------------------
C     
                silpro(ji,jj,jk)=silpro(ji,jj,jk)*trn(ji,jj,jk,jpsil)
     &              /(trn(ji,jj,jk,jpsil)+xksi2)
C     
C     Account for changes in total co2, alkalinity and [po4]
C     due to uptake by poc-/caco3 producing organismS
C     ------------------------------------------------------
C     
C     Evolution of silicates
C     ----------------------
C     
                silpro(ji,jj,jk)=amin1(silpro(ji,jj,jk),
     &              (trn(ji,jj,jk,jpsil)-0.05E-6))    
                trn(ji,jj,jk,jpsil) =trn(ji,jj,jk,jpsil)
     $              -SILPRO(ji,jj,jk)
C     
C     Consumption of Total (12C)O2
C     ----------------------------
C     
                trn(ji,jj,jk,jpdic) = trn(ji,jj,jk,jpdic)
     &              -prodt(ji,jj,jk)-prcaca(ji,jj,jk)
C     
#    if defined key_trc_biohamocc13
C     Consumption of Total (13C)O2
C     ----------------------------
                trn(ji,jj,jk,jp13c) = trn(ji,jj,jk,jp13c)-
     &              pdb*plafr13*prodt(ji,jj,jk)
     &              -pdb*prcaca(ji,jj,jk)
#    endif
C     
C     Consumption of alkalinity due to ca++ uptake and increase 
C     of alkalinity due to nitrate consumption during organic 
C     soft tissue production
C     ---------------------------------------------------------
C     
                trn(ji,jj,jk,jptal) = trn(ji,jj,jk,jptal)+
     &              rno3*prodt(ji,jj,jk)-two*prcaca(ji,jj,jk)
C     
            ENDIF
          END DO
        END DO
      ENDDO
C     
#endif
C     
      RETURN
      END



