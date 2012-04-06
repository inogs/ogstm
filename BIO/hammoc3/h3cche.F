      SUBROUTINE h3cche
#if defined key_passivetrc && defined key_trc_hamocc3
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE h3cche
CCC                     ******************
CCC
CCC     PURPOSE.
CCC     --------
CCC          *H3CCHE* SETS CHEMICAL CONSTANTS 
CCC
CCC
CC    INTERFACE.
CC     ----------
CC          *CALL* *H3CCHE*
CC
CC     METHOD.
CC     -------
CC          1) SET CONSTANTS FOR CARBONATE CHEMISTRY AS DESCRIBED IN
CC             IN BROECKER ET AL. (1982, GEOSECS) AND EDMOND A. GIESKES
CC             (1970)
CC          2) INITIATE [CO3--] AND PH-VALUE BY ITERATION
CC             (NEWTON-RAPHSON METHOD FOR SOLVING NONLINEAR SIMULTANEOUS
CC              EQUATIONS, SEE E.G. SCARBOROUGH, J. (1958))
CC
CC
CC     EXTERNALS.
CC     ----------
CC          *RHO*     - HALF PRECISION FUNCTION, EQ. OF STATE OF
CC                      SEAWATER
CC
CC     REFERENCE.
CC     ----------
CC
CC          BERNER, R. A. (1976)
CC          THE SOLUBILITY OF CALCITE AND ARAGONITE IN SEA WATER
CC          AT ATMOSPHERIC PRESSURE AND 34.5 O/OO SALINITY.
CC          AMERICAN JOURNAL OF SCIENCE, VOL. 276, 713-730.
CC          (K'SP(ARAGONITE)=1.45 K'SP(CALCITE))
CC
CC          BROECKER, W.S., D.W. SPENCER, AND H. CRAIG (1982)
CC          GEOSECS PACIFIC EXPEDITION. VOL. 3.. HYDROGRAPHIC DATA
CC          1973-1974, SUPERINTENDANT OF DOCUMENTS, U.S. GOVERNMENT
CC          PRINTING OFFICE, WASHINGTON, D.C., 137 PP..
CC
CC          CULBERSON, C.H., AND R.M. PYTKOWICZ (1968)
CC          EFFECT ON PRESSURE ON CARBONIC ACID, BORIC ACID AND THE PH
CC          IN SEA WATER.
CC          LIMNOLOGY AND OCEANOGRAPHY, VOL. 13, 403-417.
CC
CC          DICKSON, A.G., AND J.P. RILEY (1979)
CC          THE ESTIMATION OF ACID DISSOCIATION CONSTANTS IN SEAWATER
CC          MEDIA FROM POTENTIOMETRIC TITRATIONS WITH STRONG BASE.
CC          I. THE IONIC PRODUCT OF WATER - KW.
CC          MARINE CHEMISTRY, VOL. 7, 89-99.
CC
CC          EDMOND, J.M., AND J.M.T.M. GIESKES (1970)
CC          ON THE CALCULATION OF THE DEGREE OF SATURATION OF SEA WATER
CC          WITH RESPECT TO CALCIUM CARBONATE UNDER IN SITU CONDITIONS.
CC          GEOCHIM. ET COSMOCHIM. ACTA, 34, 1261-1291.
CC
CC          INGLE, S.E. (1800)
CC          SOLUBILITY OF CALCITE IN THE OCEAN.
CC          MARINE CHEMISTRY, VOL. 3, 301-319.
CC
CC          INGLE, S.E., C.H. CULBERSON, J.E. HAWLEY, AND R.M. PYTKOWICZ
CC          (1973) THE SOLUBILITY OF CALCITE IN SEAWATER AT ATMOSPHERIC
CC          PRESSURE AND 35 O/OO SALINITY.
CC          MARINE CHEMISTRY, VOL. 1, 295-307.
CC
CC          RILEY, J. P., AND G. SKIRROW, EDS. (1965)
CC          CHEMICAL OCEANOGRAPHY. VOL. 1, 712 PP., ACADEMIC PRESS,
CC          LONDON A. NEW YORK.
CC
CC          SCARBOROUGH, J. (1958) NUMERICAL MATHEMATICAL ANALYSIS.
CC          OXFORD UNIVERSITY PRESS, LONDON, 4TH ED., 576 PP..
CC
CC          WEISS, R. F. (1970) THE SOLUBILITY OF NITROGEN
CC          OXYGEN AND ARGON IN WATER AND SEAWATER.
CC          DEEP-SEA RESEARCH, VOL. 17, 721-735.
CC
CC          WEISS, R. F. (1974)
CC          CARBON DIOXIDE IN WATER AND SEAWATER: THE SOLUBILITY OF A
CC          NON IDEAL GAS. MARINE CHEMISTRY, VOL. 2, 203-215.
CC
CC          WOOSTER, W.S., A.J. LEE, AND G. DIETRICH (1969)
CC          REDEFINITION OF SALINITY. Z. GEOPHYS., VOL.35, 611-613.
CC
CC          BROECKER, W.S., D.W. SPENCER, AND H. CRAIG (1982)
CC          GEOSECS PACIFIC EXPEDITION. VOL. 3.. HYDROGRAPHIC DATA
CC          1973-1974, SUPERINTENDANT OF DOCUMENTS, U.S. GOVERNMENT
CC          PRINTING OFFICE, WASHINGTON, D.C., 137 PP..
CC
CC
CC*     VARIABLE           TYPE    PURPOSE.
CC      --------           ----    --------
CC
CC      *PRES              REAL(8)    APPROXIMATE PRESSURE AT DEPTH OF U-POINTS
CC                                 IN BAR, DUMMY VARIABLE
CC      *TC*               REAL(8)    TEMPERATURE AT OCEAN GRID POINTS (DEG C),
CC                                 DUMMY VARIABLE
CC      *CL*               REAL(8)    CHLORINITY (CL(O/OO)=S(O/OO)/1.80655)
CC                                 AFTER WOOSTER ET AL., 1969
CC                                 (C.F. KALLE/DIETRICH , P. 60)
CC      *AKW*              REAL(8)    KW, H2O DISSOC. CONSTANT, LIT ?
CC      *H*                REAL(8)    [H+], DUMMY VARIABLE
CC      *rrr               REAL(8)    [CO3--] [MOLE/L], DUMMY VARIABLE
CC      *C*                REAL(8)    GIVEN [SUM(12C)O2] [MOLE/L], DUMMY VARIABLE
CC      *A*                REAL(8)    ALKALINITY [EQV/L] AS FUNCTION OF [CO3--]
CC                                 AND [H+], DUMMY VARIABLE
CC   MODIFICATIONS:
CC   --------------
CC      original :      1988 E. Maier-Reimer
CC      additions :     1998 O. Aumont
CC      modifications : 1999 C. Le Quere
CC----------------------------------------------------------------------
CC parameters and commons
CC ======================

       USE myalloc
        IMPLICIT NONE

CC----------------------------------------------------------------------
CC local declarations
CC ==================
C
      INTEGER ji, jj, jk
      REAL(8) tkel, sal
      REAL(8) pres, tc, r, cl
      REAL(8) akw, a, c, h, akb
      REAL(8) zsqrt, ztr, zlogt
      REAL(8) zqtt,qtt2,sal2,sal15
C
C* 21. CHEMICAL CONSTANTS - SURFACE LAYER
C ---------------------------------------
C
C vertical slab
C =============
C
      DO jj = 1,jpj
C
        DO ji = 1,jpi
C
C* 21.1 SET ABSOLUTE TEMPERATURE
C ------------------------------
C
          tkel = tn(ji,jj,1)+temzer
          qtt = tkel*perc
          qtt2=qtt**2
          sal = sn(ji,jj,1) + (1.-tmask(ji,jj,1))*35.
          zqtt=log(qtt)

C
C* 21.2 CHLORINITY (WOOSTER ET AL., 1969)
C ---------------------------------------
C
          cl = sal*salchl
C
C* 21.3 LN(K0) OF SOLUBILITY OF CO2 (EQ. 12, WEISS, 1974)
C -------------------------------------------------------
C
          cek0 = c00+c01/qtt+c02*zqtt+sal*(c03+c04*qtt+c05*qtt2)
C
C* 21.6 LN(K0) OF SOLUBILITY OF O2 and N2 (EQ. 4, WEISS, 1970)
C ------------------------------------------------------------
C
          oxy = ox0+ox1/qtt+ox2*zqtt+sal*(ox3+ox4*qtt+ox5*qtt2)
C
C* 21.9 CO2 SOLUBILITY IN SEAWATER (WEISS, 1974, CF. EQ. 12)
C ----------------------------------------------------------
C
          ak0 = exp(cek0)*smicr
C
C* 21.10 DENSITY OF SEAWATER AND TOTAL BORATE IN MOLE/
C ----------------------------------------------------
C
          rrr = rhopn(ji,jj,1) *thousi
          bor = bor1*rrr*cl*bor2
C
C* 21.11 SET CHEMICAL CHEMICAL CONSTANTS
C --------------------------------------
C
          chemc(ji,jj,1) = ak0
          chemc(ji,jj,2) = bor
C
C* 21.12 O2 SOLUBILITY IN SEAWATER (WEISS, 1970, CF. EQ. 4)
C ---------------------------------------------------------
C
          chemc(ji,jj,3) = exp(oxy)*oxyco
C
        ENDDO
      END DO
C
C*  22. CHEMICAL CONSTANTS - DEEP OCEAN
C -------------------------------------
C
      DO jk = 1,jpk
C
C* 22.1 APPROX. SEAWATER PRESSURE AT U-POINT DEPTH (BAR)
C ------------------------------------------------------
C
        pres = 1.025e-1*gdept(jk)
C
        DO jj=1,jpj
          DO ji = 1,jpi
C
C* 22.2 SET LIMITS FOR SEAWATER TEMP. AND SALINITY
C              (THIS IS DONE TO AVOID COMPUTATIONAL CRASH AT DRY
C               POINTS DURING CALCULATION OF CHEMICAL CONSTANTS)
C
C* 22.3 SET [H+] (FIRST GUESS)
C ----------------------------
C
            hi(ji,jj,jk) = tenm7
C
C* 22.4 SET ABSOLUTE TEMPERATURE
C ------------------------------
C
            tkel   = tn(ji,jj,jk)+temzer
            qtt    = tkel*perc
            sal    = sn(ji,jj,jk) + (1.-tmask(ji,jj,jk))*35.
            zsqrt  = sqrt(sal)
            sal15 = sal**1.5
            zlogt  = log(tkel)
            ztr    = 1./tkel
C
C* 22.5 CHLORINITY (WOOSTER ET AL., 1969)
C ---------------------------------------
C
            cl = sal*salchl
C
C* 22.6 LN(K0) OF SOLUBILITY OF CO2 (EQ. 12, WEISS, 1974)
C -------------------------------------------------------
C
            cek0 = c00+c01/qtt+c02*alog(qtt)+
     &          sal*(c03+c04*qtt+c05*qtt**2)
C
C  COEFFICIENT OCMIP 
C ------------------
C
            ckb = (cb0+cb1*zsqrt+cb2*sal+cb3*sal15+cb4*sal**2)*ztr
     $          +(cb5+cb6*zsqrt+cb7*sal)+
     $          (cb8+cb9*zsqrt+cb10*sal)*zlogt+cb11*zsqrt*tkel
            ck1 = c10*ztr+c11+c12*zlogt+(c13*ztr+c14)*zsqrt+
     $          c15*sal+c16*sal15+log(1.+c17*sal)
            ck2 = c20*ztr+c21+c22*zlogt+(c23*ztr+c24)*zsqrt+c25*sal
     $          +c26*sal15+log(1.+c27*sal)
C
C* 21.5 PKW (H2O) (DICKSON AND RILEY, 1979)
C -----------------------------------------
C
            ckw = cw0*ztr+cw1+cw2*zlogt+(cw3*ztr+cw4+cw5*zlogt)*
     $          zsqrt+cw6*sal
C
C* 22.8 LN(K0) OF SOLUBILITY OF O2 (EQ. 4, WEISS, 1970)
C -----------------------------------------------------
C
            oxy = ox0+ox1/qtt+ox2*alog(qtt)+sal*(ox3+ox4*qtt+ox5*qtt**2)
C
C* 22.9 K1, K2 OF CARBONIC ACID, KB OF BORIC ACID, KW (H2O) (LIT.?)
C -----------------------------------------------------------------
C
            ak1 = exp(ck1)
            ak2 = exp(ck2)
            akb = exp(ckb)
            akw3(ji,jj,jk) = exp(ckw)
C
C*22.10 APPARENT SOLUBILITY PRODUCT K'SP OF CALCITE IN SEAWATER
C       (S=27-43, T=2-25 DEG C) AT pres =0 (ATMOSPH. PRESSURE)
C       (INGLE, 1800, EQ. 6)
C -------------------------------------------------------------
C
            aksp0 = 1.E-7*(akcc1+akcc2*sal**(third)+akcc3*alog10(sal)
     &          +akcc4*tkel**2)
C
C* 22.11 FORMULA FOR CPEXP AFTER EDMOND AND GIESKES (1970)
C
C        (REFERENCE TO CULBERSON AND PYTKOQICZ (1968) AS MADE
C        IN BROECKER ET AL. (1982) IS INCORRECT; HERE RGAS IS
C        TAKEN TENFOLD TO CORRECT FOR THE NOTATION OF pres  IN
C        DBAR INSTEAD OF BAR AND THE EXPRESSION FOR CPEXP IS
C        MULTIPLIED BY LN(10.) TO ALLOW USE OF EXP-FUNCTION
C        WITH BASIS E IN THE FORMULA FOR AKSPP (CF. EDMOND
C        AND GIESKES (1970), P. 1285 AND P. 1286 (THE SMALL
C        FORMULA ON P. 1286 IS RIGHT AND CONSISTENT WITH THE
C        SIGN IN PARTIAL MOLAR VOLUME CHANGE AS SHOWN ON
C        P. 1285))
C -----------------------------------------------------------
C
            cpexp = pres /(rgas*tkel)
C
C* 22.12 KB OF BORIC ACID, K1,K2 OF CARBONIC ACID PRESSURE
C        CORRECTION AFTER CULBERSON AND PYTKOWICZ (1968)
C        (CF. BROECKER ET AL., 1982)
C --------------------------------------------------------
C
            tc = tn(ji,jj,jk) + (1.-tmask(ji,jj,jk))*20.
            akb3(ji,jj,jk) = akb*exp(cpexp*(devkb-devkbt*tc))
            ak13(ji,jj,jk) = ak1*exp(cpexp*(devk1-devk1t*tc))
            ak23(ji,jj,jk) = ak2*exp(cpexp*(devk2-devk2t*tc))
C
C  22.13 APPARENT SOLUBILITY PRODUCT K'SP OF CALCITE (OR ARAGONITE)
C        AS FUNCTION OF PRESSURE FOLLWING EDMOND AND GIESKES (1970)
C        (P. 1285) AND BERNER (1976)
C -----------------------------------------------------------------
C
            aksp(ji,jj,jk) = aracal*aksp0*exp(cpexp*(devks-devkst*tc))
C
C* 22.14 DENSITY OF SEAWATER AND TOTAL BORATE CONCENTR. [MOLES/L]
C ---------------------------------------------------------------
C
            rrr = rhopn(ji,jj,jk)*thousi
            bor = bor1*rrr*cl*bor2
            borat(ji,jj,jk) = bor
C
          ENDDO
        ENDDO
      END DO
C      
#endif
C
      RETURN
      END

