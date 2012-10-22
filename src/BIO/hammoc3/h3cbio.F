      SUBROUTINE h3cbio
#if defined key_passivetrc && defined key_trc_hamocc3
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE h3cbio
CCC                     ******************
CCC
CCC
CCC     PURPOSE.
CCC     --------
CCC          *H3CBIO* MODELS PRODUCTION OF BIOGENIC MATTER (POC ''SOFT
CCC                   TISSUE'' AND CACO3 PARTICLES ''HARD PARTS'')
CCC                   AND ITS DISTRIBUTION IN WATER COLUMN
CCC
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
CC          GLOBAL TRANSPORT", KANDERSTEG, 2 TO 6 SEPTEMBER 1985,
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
CC      original      : 1988-07 E. MAIER-REIMER      MPI HAMBURG 
CC      additions     : 1998    O. Aumont
CC      modifications : 1999    C. Le Quere
CC      modifications : 2001    O. Aumont
CC ---------------------------------------------------------------------------
CC parameters and commons
CC ======================
       USE myalloc
        IMPLICIT NONE

CC----------------------------------------------------------------------
CC local declarations
CC ==================
C
      INTEGER ji,jj,jk,ktask
      REAL(8) propro, brupro, silpot, expofa
      REAL(8) calred, calpot, topot, silfra
C
C ============================================================
C* 2. CALCULATE PRODUCTION OF SOFT TISSUE AND RELATED CACO3
C     HARD PARTS USING MICHAELIS MENTON KINETICS (CF. DUGDALE,
C     1967; PARSONS AND TAKAHASHI, 1973)
C ============================================================
C
C* 1. SET HALF PRECISION CONSTANTS
C -----------------------------------------------------------------
C
      zero   = 0.
      one    = 1.
      two    = 2.
C
      CALL h3copt
C
C vertical slab
C ===============
C
C
        DO jk = 1,jpkb
          DO jj = 1,jpj
            DO ji = 1,jpi
C
C* 2.2  PRIMARY PRODUCTION OF POC (SOFT TISSUE) USING
C       MICHAELIS-MENTON NUTRIENT UPTAKE KINETICS [MOLE/L]
C --------------------------------------------------------
C
          propro = strn(ji,jj)*0.25
          IF (propro.le.0.) propro = 0.
C
C  for T lt -2. no production 
C ---------------------------
C
          IF (tn(ji,jj,jk).lt.-2.) THEN
            prlat(ji,jj) = 0.
          ELSE
            prlat(ji,jj) = propro*((tn(ji,jj,jk)+2.)/(tn(ji,jj,jk)+10.))
     $      *tmask(ji,jj,jk)*(1.-freeze(ji,jj))
C
#    if defined key_off_degrad
     $      *facvol(ji,jj,jk)
#    endif
          ENDIF
C
          IF (hmld(ji,jj).gt.50.) THEN
             prlat(ji,jj) = prlat(ji,jj)*(50./hmld(ji,jj))              
          ENDIF
C
C* 2.3  POC PRODUCTION RELATED BIOGENIC CACO3 PARTICLE
C       PRODUCTION [MOLE/L]
C ----------------------------------------------------
C
          prorca(ji,jj,jk) = prlat(ji,jj)*trn(ji,jj,jk,jppo4)**2
     &              /(conc0+trn(ji,jj,jk,jppo4))*rfact

          brupro = grosip*prorca(ji,jj,jk)
          tra(ji,jj,jk,jppo4) = tra(ji,jj,jk,jppo4)-
     &                                    prorca(ji,jj,jk)*rfactr
          tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy)+
     &                                    o2ut*prorca(ji,jj,jk)*rfactr
          silpot = min(brupro,trn(ji,jj,jk,jpsil))
          expofa = exp(0.1*(tn(ji,jj,jk)-10.))
          calred = expofa/(1.+expofa)
          calpot = calred*prorca(ji,jj,jk)
          topot  = silpot+calpot+1.E-20
          silfra = silpot/topot
C ENDE ERSATZ
          prcaca(ji,jj,jk) = caco3r*calpot*silic0/(silic0+silpot) 
          silpro(ji,jj,jk) = silpot*silfra
     $              *trn(ji,jj,jk,jpsil)/(trn(ji,jj,jk,jpsil)+silic0)
C
C* 5.1  SURFACE LAYER C13/C12, C14/C12 TOTAL CO2 RATIOS
C -----------------------------------------------------
C
            tra(ji,jj,jk,jpsil) = tra(ji,jj,jk,jpsil)-
     &                                     silpro(ji,jj,jk)*rfactr
C
C* 5.2  CONSUMPTION OF TOTAL (12C)O2
C ----------------------------------
C
            tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) -
     &                (prorca(ji,jj,jk)+prcaca(ji,jj,jk)) * rfactr
C
C* 5.3  CONSUMPTION OF ALKALINITY DUE TO CA++ UPTAKE
C       AND INCREASE OF ALKALINITY DUE TO NITRATE
C       CONSUMPTION DURING ORGANIC SOFT TISSUE PRODUCTION
C -------------------------------------------------------
C
            tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal)+
     &           (RNO3*PRORCA(ji,jj,jk)-TWO*PRCACA(ji,jj,jk))*rfactr
C
C* 5.4  CONSUMPTION OF PHOSPHATE
C       PHOSPHATE REDUCTION HAPPENS IN SUBROUTINE ADVECT
c         trn(ji,jj,jk,jppo4) = trn(ji,jj,jk,jppo4)-PRORCA(ji,jj,jk)
C ------------------------------------------------------
C
#    if defined key_trc_biohamocc13
C
C* 5.6  CONSUMPTION OF TOTAL (13C)O2
C ----------------------------------
C
            tra(ji,jj,jk,jp13c) = tra(ji,jj,jk,jp13c) -
     &       (pdb*plafr13*prorca(ji,jj,jk)+pdb*prcaca(ji,jj,jk))*rfactr
C
#    endif
C
          ENDDO
        END DO
      END DO
C
#     if defined key_trc_diaadd
           do jj=1,jpj
             do ji=1,jpi
              trc2d(ji,jj,8)=tn(ji,jj,1)
              trc2d(ji,jj,9)=hmld(ji,jj)
             end do
           end do

#     endif
#endif
      RETURN
      END

