      SUBROUTINE h3cflx
#if defined key_passivetrc && defined key_trc_hamocc3
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE h3cflx
CCC                     ******************
CCC
CCC
CC     PURPOSE.
CC     --------
CC          *H3CFLX* CALCULATES GAS EXCHANGE AND CHEMISTRY AT SEA SURFACE
CC
CC     METHOD.
CC     -------
CC          SOLVING SYSTEM OF TWO NON-LINEAR SIMULTANEOUS EQUATIONS
CC     FOR [H2CO3] AND [H+] 
CC
CC     EXTERNALS.
CC     ----------
CC          NONE.
CC
CC     REFERENCE.
CC     ----------
CC
CC          BROECKER, W.S., AND T.-H. PENG (1982)
CC          TRACERS IN THE SEA.
CC          ELDIGIO PRESS, LAMONT-DOHERTY GEOLOGICAL OBSERVATORY,
CC          PALISADES, N.Y., 690 PP..
CC
CC          MOOK, W.G. (1986)
CC          13C IN ATMOSPHERIC CO2.
CC          NETHERLANDS JOURNAL OF SEA RESEARCH, 20(2/3): 211-223.
CC
CC          SCARBOROUGH, J. (1958) NUMERICAL MATHEMATICAL ANALYSIS.
CC          OXFORD UNIVERSITY PRESS, LONDON, 4TH ED., 576 PP..
CC
CC*     VARIABLE           TYPE    PURPOSE.
CC      --------           ----    --------
CC
CC      *ZEROAB*           REAL(8)    273.15 (0 DEG C IN DEG K), HALF PREC.
CC      *ONR*              REAL(8)    1., HALF PRECISION
CC      *TWR*              REAL(8)    2., HALF PRECISION
CC      *FOURR*            REAL(8)    4., HALF PRECISION
CC      *EVFR00*           REAL(8)    COEFFICIENT OF TEMPERATURE DEPENDENT
CC                                 ISOTOPIC CARBON FRACTIONATION FACTOR
CC                                 FOR EVAPORATION (MOOK, 1986)
CC      *EVFR01*           REAL(8)    COEFFICIENT OF TEMPERATURE DEPENDENT
CC                                 ISOTOPIC CARBON FRACTIONATION FACTOR
CC                                 FOR EVAPORATION (MOOK, 1986)
CC      *XN*               REAL(8)    SURFACE [H+] (EXPRESSED AS
CC                                 SQRT(AK1*AK2)/[H+])
CC                                 AT END OF ONE ITERATION STEP
CC      *CTO*              REAL(8)    ACTUAL [SUM(C-12)O2], DUMMY VARIABLE
CC      *CT1*              REAL(8)    ACTUAL [SUM(C-12)O2], DUMMY VARIABLE
CC      *ALKA*             REAL(8)    ACTUAL [ALK] [EQV/L], DUMMY VARIABLE
CC      *AKB*              REAL(8)    1. DISSOC. CONSTANT OF BORIC ACID
CC      *H2CO3(jpi,jpj)*     REAL(8)    SURFACE [H2CO3] [MOLE/L]
CC      *TWR*              REAL(8)    2., HALF PRECISIONION
CC      *FOURR*            REAL(8)    4., HALF PRECISIONION
CC      *YEAR*             REAL(8)    NUMBER OF SECONDS IN 1 YEAR
CC      *KRORR*         INTEGER    COUNTS ITERATIONS IN NEWTON-RAPHSON ALGO-
CC                                 RITHM
CC      *alpco2*            REAL(8)    AK0, SOLUBILITY OF CO2 IN SEAWATER, DUMMY
CC                                 VARIABLE
CC      *BT*               REAL(8)    TOTAL BORATE [MOLES/L]
CC                                 ([B(OH)3]+[B(OH)4-]), DUMMY VARIABLE
CC      *CT1*              REAL(8)    ACTUAL VALUE OF INORG. C, DUMMY VARIABLE
CC      *X1*               REAL(8)    [H+] EXPRESSED AS SQRT(AK1*AK2)/[H+],
CC                                 DUMMY VARIABLE
CC      *PCO2*             REAL(8)    PCO2 OF SURFACE WATER [PPM], DUMMY VARIABLE
CC      *PA*               REAL(8)    ATMOSPHERIC PCO2 [PPM], DUMMY VARIABLE
CC      *FUGAC*            REAL(8)    ENHANCEMENT FACTOR FOR GAS EXCHANGE FLUXES
CC      *FLD*              REAL(8)    [CO2] FLUX ATMOSPHERE -> OCEAN
CC                                 IN [PPM/TIME STEP]
CC      *FLU*              REAL(8)    [CO2] FLUX OCEAN -> ATMOSPHER
CC                                 IN [PPM/TIME STEP]
CC      *RELW13*           REAL(8)    RATIO [SUM((13C)O2)]/[SUM((12C)O2)]
CC                                 IN SURFACE LAYER
CC      *RELA13*           REAL(8)    RATIO [SUM((13C)O2)]/[SUM((12C)O2)]
CC                                 IN THE ATMOSPHERE
CC      *RELW14*           REAL(8)    RATIO [SUM((14C)O2)]/[SUM((12C)O2)]
CC                                 IN SURFACE LAYER
CC      *RELA14*           REAL(8)    RATIO [SUM((14C)O2)]/[SUM((12C)O2)]
CC                                 IN THE ATMOSPHERE
CC
CC   MODIFICATIONS:
CC   --------------
CC      original      : 1988-07 E. MAIER-REIMER      MPI HAMBURG
CC      additions     : 1998    O. Aumont
CC      modifications : 1999    C. Le Quere
CC     -----------------------------------------------------------------
CC  parameters and commons
CC ======================

       USE myalloc
       IMPLICIT NONE

CC----------------------------------------------------------------------
CC local declarations
CC ==================
C
      INTEGER ji, jj, it
      REAL(8) zexp1, zexp2
      REAL(8) a1, a2, a3, b2, b3, ttc, ws, krorr, alpco2
      REAL(8) pa, fld, flu, conveg
      REAL(8) oxy16, voro16, flu16
      REAL(8) relw13, rela13, fra13, frw13
      REAL(8) h,ah2,bot
      REAL(8) ct1,a, alka
      REAL(8) schmitt, caralk, bicarb
C
C  1. ASSIGNATION TO EXPONENTS IN THE LISS AND MERLIVAT
C     FORMULATION OF THE GAS EXCHANGE RATE
c -----------------------------------------------------
C
      zexp1 = -2./3.
      zexp2 = -1./2.  
      a1    = 0.17
      a2    = 2.85
      a3    = 5.90
      b2    = 9.65
      b3    = 49.3      
C
C  2.1 INITIALISATION OF [H+]
C  --------------------------
C
      DO jj=1,jpj
        DO ji=1,jpi
          caralk = trn(ji,jj,1,jptal)-
     &                 borat(ji,jj,1)/(1.+1E-8/akb3(ji,jj,1))
          co3(ji,jj,1) = caralk-trn(ji,jj,1,jpdic)
     &        +(1.-tmask(ji,jj,1))*.5e-3
          bicarb = (2.*trn(ji,jj,1,jpdic)-caralk)
          hi(ji,jj,1) = ak23(ji,jj,1)*bicarb/co3(ji,jj,1)
        END DO 
      END DO 
      
C
C* 2.2 SURFACE CHEMISTRY (PCO2 AND [H+] IN
C     SURFACE LAYER); THE RESULT OF THIS CALCULATION
C     IS USED TO COMPUTE AIR-SEA FLUX OF CO2
C ---------------------------------------------------
C
      DO krorr = 1,30
C
        DO jj = 1,jpj
          DO ji = 1,jpi
C
C* 2.3 TOTAL BORATE
C -----------------
C
            bot = chemc(ji,jj,2)
            ct1  = trn(ji,jj,1,jpdic)
            alka = trn(ji,jj,1,jptal) 
            h = hi(ji,jj,1)+(1.-tmask(ji,jj,1))*1.e-9
            h = amax1(hi(ji,jj,1),1.E-10)
            
C
C* 2.4 CALCULATE [ALK]([CO3--], [HCO3-])
C ------------------------------------
C
            a=alka-
     &          (akw3(ji,jj,1)/h-h+bot/(1.+h/akb3(ji,jj,1)))
C
C* 2.5 CALCULATE [H+] AND [H2CO3]
C -----------------------------------------
C
            ah2=sqrt((ct1-a)**2+4*(a*ak23(ji,jj,1)/ak13(ji,jj,1))
     &          *(2*ct1-a)) 
            ah2=0.5*ak13(ji,jj,1)/a*((ct1-a)+ah2)
            hi(ji,jj,1)  = ah2
            h2co3(ji,jj) = (2*ct1-a)/(2.+ak13(ji,jj,1)/ah2)
            conveg=((ah2-hi(ji,jj,1))/hi(ji,jj,1))**2
     $          *tmask(ji,jj,1)
            rconvs = rconvs+conveg
          END DO
        END DO
C
C  2.6 TEST CONVERGENCE
C ---------------------
C
        IF (rconvs.LE.1.E-7) EXIT
C
      END DO
C
C
C  2.7 CHECK ITERATION
C  --------------------
C
C     IF (krorr.GE.30) THEN 
C         WRITE(numout,*) 'in h3cflx h2co3 have not converged '
C         CALL FLUSH_(numout)
C         STOP 'h3cflx'
C     END IF 
C
C
C 3. COMPUTE FLUXES
C --------------
      
C
C 3.1 FIRST COMPUTE GAS EXCHANGE COEFFICIENTS
C -------------------------------------------
C
      
      DO jj = 1,jpj
        DO ji = 1,jpi
C
          ttc = min(39.,tn(ji,jj,1))
          schmitt= 2073.1-125.62*ttc+3.6276*ttc**2-0.043126*ttc**3
C
C 3.2 COMPUTE GAS EXCHANGE FOR CO2
C --------------------------------
C
          if (igaswind.eq.1) then
              ws = vatm(ji,jj)
C
C 3.3 this is Wanninkopf(1992) equation 8 (with chemical enhancement), 
C     in cm/h
C --------------------------------------------------------------------
C
              kgwanin(ji,jj) = (0.3*ws*ws + 2.5*(0.5246+ttc*(0.016256+
     &            ttc*0.00049946)))*sqrt(660./schmitt)
C
C 3.4 CONVERT TO M/S
C ------------------
C
              kgwanin(ji,jj) = kgwanin(ji,jj)/100./3600.
C
C 3.5 convert to mol/m2/s/uatm, alpco2(chemc(ji,jj,1)) is in 
C      mol/L/uatm and apply ice cover
C -----------------------------------------------------------
C
              kgwanin(ji,jj) = kgwanin(ji,jj)*chemc(ji,jj,1)*1.e3 * 
     &            (1-freeze(ji,jj))
C
C 3.6 CLIMATOLOGICAL CASE: GAS EXCHANGE IS NOT COMPUTED BUT
C     READ IN A FILE
C ---------------------------------------------------------
C
           ELSE IF (igaswind.EQ.2) THEN
C
C 3.7 CONVERT TO M/S
C ------------------
C
               kgwanin(ji,jj) = kgwanin(ji,jj)/raass
     &             *(1-freeze(ji,jj))
C
C
C 3.8 COUPLED RUN CASE. WIND FIELD IS COMPUTED FROM WIND STRESS
C --------------------------------------------------------------
C
           ELSE IF (igaswind.EQ.3) THEN
               
C 3.9 LB zonal stress = 0 on east side of bassins
C     LB uses stress from coupled runs to compute kgwanin
C --------------------------------------------------------
C
               ws=sqrt((sqrt((taux(ji,jj)**2)+(tauy(ji,jj))**2))
     $             /(1.25e-3))
               kgwanin(ji,jj)=(0.3*(sqrt((taux(ji,jj)**2)+
     $             (tauy(ji,jj)) **2))/(1.25e-3)  + 2.5*(0.5246 +
     $             ttc*(0.016256 +ttc*0.00049946)))*sqrt(660./schmitt)
C
C 3.10 CONVERT TO M/S
C --------------------
C
               kgwanin(ji,jj) = kgwanin(ji,jj)/360000.
C
C 3.11 convert to mol/m2/s/uatm, alpco2(chemc(ji,jj,1)) is in mol/L/uatm
C      and apply ice cover
C ----------------------------------------------------------------------
C
               kgwanin(ji,jj)=kgwanin(ji,jj)*chemc(ji,jj,1)*1.e3
     &             *(1-freeze(ji,jj))
           ENDIF
         END DO
       END DO
C
C
C 3.12 COMPUTE GAS EXCHANGE COEFFICIENT FO O2 FROM LISS AND 
C      MERLIVAT EQUATIONS
C ---------------------------------------------------------
C
       DO jj = 1,jpj
         DO ji=1,jpi
C
           IF (igaswind.EQ.3) then
               ws=sqrt((sqrt((taux(ji,jj)**2)+(tauy(ji,jj))**2))
     $             /(1.25e-3))
           ELSE
               ws = vatm(ji,jj)
           ENDIF
           
           ttc = min(39.,tn(ji,jj,1))
           SChmitt = 1953.4-128.0*ttc+3.9918*ttc**2-0.050091*ttc**3
C
C 
C
           IF (ws.LE.3.6) THEN
               fugaci(ji,jj) = (a1*ws)*(schmitt/660.)**zexp1
           ENDIF
           IF ((ws.GT.3.6).AND.(ws.LE.13.)) THEN
               fugaci(ji,jj) = (a2*ws-b2)*(schmitt/660.)**zexp2
           ENDIF
           IF (ws.GT.13.) THEN
               fugaci(ji,jj) = (a3*ws-b3)*(schmitt/660.)**zexp2
           ENDIF
C
C convert to cm and apply ice cover
C
           fugaci(ji,jj) = fugaci(ji,jj)/100./(3600.)*
     $         (1-freeze(ji,jj))*tmask(ji,jj,1)
C
#    if defined key_off_degrad
           fugaci(ji,jj) = exp(-rfact*fugaci(ji,jj)
     $         *facvol(ji,jj,1)/e3t(1))
#    else
           fugaci(ji,jj) = exp(-rfact*fugaci(ji,jj)
     $         /e3t(1))
#    endif
           
         ENDDO
       ENDDO
C
       DO jj = 1,jpj
         DO ji = 1,jpi
C
C Fix atmospheric CO2
C -------------------
C
C for climate runs
           pa = atcco2 
C
C Compute CO2 flux for the sea and air
C ------------------------------------
C
           alpco2 = chemc(ji,jj,1)
           fld = pa*tmask(ji,jj,1)*kgwanin(ji,jj)
           flu = h2co3(ji,jj)/alpco2*tmask(ji,jj,1)*kgwanin(ji,jj)
           tra(ji,jj,1,jpdic)= tra(ji,jj,1,jpdic)+(fld-flu)/1000./e3t(1)
CMAFOA faire qcumul sur domaine exact
           qcumul(1)=qcumul(1)+(fld-flu)*rfact*e1t(ji,jj)
     $         *e2t(ji,jj)*tmask(ji,jj,1)
C
#    if defined key_trc_biohamocc13
C
C Compute C13 flux
C ----------------
C
          relw13 = trn(ji,jj,1,jp13c)/
     &       (trn(ji,jj,1,jpdic)+(1.-tmask(ji,jj,1))*1.e-20)
          rela13 = pdb*(1.-6.4/1000.)
          fra13  = 1.+(-10.66+0.1196*tn(ji,jj,1)-0.0003095*
     &          tn(ji,jj,1)**2)/1000.
          frw13  = 1.
          tra(ji,jj,1,jp13c) = tra(ji,jj,1,jp13c)+
     &       (fld*frw13*rela13-flu*relw13*fra13)/1000./e3t(1)
#    endif
C
C Compute O2 flux 
C ---------------
C
          oxy16 = trn(ji,jj,1,jpoxy)
          voro16 = atcox
          flu16 = (-fugaci(ji,jj)+1)*e3t(1)
     &            *(voro16*chemc(ji,jj,3)-oxy16)*
     &            tmask(ji,jj,1)/rfact 
          tra(ji,jj,1,jpoxy) = tra(ji,jj,1,jpoxy)+flu16/e3t(1)
C
C
C
C Save diagnostics
C ----------------
C
#    if defined key_trc_diaadd
          trc2d(ji,jj,1) = (fld-flu)
          trc2d(ji,jj,2) = flu16*1000.
          trc2d(ji,jj,3) = kgwanin(ji,jj)
          trc2d(ji,jj,4) = (fld-flu)/(kgwanin(ji,jj)+1.E-15)
C
#        if defined key_trc_biohamocc13
          trc2d(ji,jj,10) = (fld*frw13*rela13-flu*relw13*fra13)
#        endif
C
#    endif
C
        END DO
      END DO
C
#endif
      RETURN
      END

