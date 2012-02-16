
       SUBROUTINE h3clys
#if defined key_passivetrc && defined key_trc_hamocc3
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE h3clys
CCC                     ******************
CCC
CCC
CCC     PURPOSE.
CCC     --------
CCC          *H3CLYS*  CALCULATES DEGREE OF CACO3 SATURATION IN THE WATER
CCC                    COLUMN, DISSOLUTION/PRECIPITATION OF CACO3 AND LOSS
CCC                    OF CACO3 TO THE CACO3 SEDIMENT POOL.
CCC
CC
CC     METHOD.
CC     -------
CC          [H+] AND [CO3--] FOR THE ACTUAL TIME STEP ARE CALCULATED
CC     BY NEWTON-RAPHSON ITERATION (E.G. SCARBOROUGH, 1958).
CC
CC     EXTERNALS.
CC     ----------
CC          NONE.
CC
CC     REFERENCE.
CC     ----------
CC
CC          SCARBOROUGH, J. (1958) NUMERICAL MATHEMATICAL ANALYSIS.
CC          OXFORD UNIVERSITY PRESS, LONDON, 4TH ED., 576 PP..
CC
CC*     VARIABLE           TYPE    PURPOSE.
CC      --------           ----    --------
CC 
CC      *NYEAR*            INTEGER COUNTS TIMESTEPS (YEARS) OF INTEGRATION
CC                           (INTEGER, INPUT)
CC      *CONVEG*           REAL(8)    CHECK FOR CONVERGENCE OF NEWTON-RAPHSON
CC                                 METHOD
CC      *KITTER*           INTEGER SETS UPPER LIMIT FOR NUMBER OF ITERATIONS
CC                                 TO DETERMINE [CO3--] AND [H+]
CC      *AKW*              REAL(8)    APPROXIMATE VALUE OF IONIC PRODUCT OF
CC                                 WATER
CC      *H*                REAL(8)    [H+], DUMMY VARIABLE
CC      *R*                REAL(8)    [CO3--] [MOLE/L], DUMMY VARIABLE
CC      *ALKA*             REAL(8)    GIVEN ALKALINITY [EQV/L], DUMMY VARIABLE
CC      *C*                REAL(8)    GIVEN [SUM(12C)O2] [MOLE/L], DUMMY VARIABLE
CC      *A*                REAL(8)    ALKALINITY [EQV/L] AS FUNCTION OF [CO3--]
CC                                 AND [H+], DUMMY VARIABLE
CC      *DELCO3*           REAL(8)    DEVIATION OF ACTUAL CACO3 CONCENTRATION FROM
CC                                 SATURATION VALUE
CC      *UNDSAT*           REAL(8)    UNDERSATURATION OF CACO3 (E.G. 3.=THREEFOLD)
CC      *EXCESS*           REAL(8)    EXCESS OF CACO3 (E.G. 3.=THREEFOLD)
CC      *DISPOT*           REAL(8)    FRACTION CACO3 (12C) THAT IS DISSOLVED
CC      *EXCE14*           REAL(8)    SUPERSATURATION IN CA(14C)O3 (E.G. 3.=
CC                                 THREEFOLD)
CC      *DISP14*           REAL(8)    FRACTION CACO3 (14C) THAT IS DISSOLVED
CC      *EXCE13*           REAL(8)    SUPERSATURATION IN CA(13C)O3 (E.G. 3.=
CC                                 THREEFOLD)
CC      *DISP13*           REAL(8)    FRACTION CACO3 (13C) THAT IS DISSOLVED
CC      *SEDLOS*           REAL(8)    FRACTION OF CACO3 IN THE BOTTOM WATER LAYER
CC                                 LOST TO THE CACO3 SEDIMENT POOL
CC      *SEDLOI*           REAL(8)    FRACTION OF CACO3 IN THE BOTTOM WATER LAYER
CC                                 WHICH REMAINS IN THE WATER COLUMN
CC
CC   MODIFICATIONS:
CC   --------------
CC      original      : 1988-07 E. MAIER-REIMER      MPI HAMBURG
CC      additions     : 1998    O. Aumont
CC      modifications : 1999    C. Le Quere
CC ---------------------------------------------------------------------------
CC parameters and commons
CC ======================

c#include "parameter.h"
c#include "common.h"
       USE myalloc

CC+CC  Implicit typing is never allowed
        IMPLICIT NONE
CC+CC  Implicit typing is never allowed

CC----------------------------------------------------------------------
CC local declarations
CC ==================
C
      INTEGER ji, jj, jk, jn
      INTEGER kitter
      REAL(8) bot, alka
      REAL(8) r, a, c
      REAL(8) delco3, excess, dispot
      REAL(8) h,remco3,ah2
      REAL(8) conveg, bicarb, caralk
C
C ------------------------------------------------------------------
C
C* 1. SET HALF PRECISION CONSTANTS
C --------------------------------
C
      zero = 0.
C
C ===========================================================
C* 2. ITERATION TO DETERMINE [CO3--] AND [H+]
C     (NEWTON-RAPHSON METHOD:
C     THE VALUES OF [SUM(CO2)] AND [ALK] ARE GIVEN,
C     DESIRED ROOTS OF [CO3--] AND [H+] FOR THAT PAIR
C     ARE DETERMINED BY SOLVING NUMERICALLY THE SYSTEM
C     OF THE TWO NONLINEAR EQUATIONS
C     1) [ALK]GIVEN      - [ALK]([CO3--],[H+])      = 0 (=F)
C     2) [SUM(CO2)]GIVEN - [SUM(CO2)]([CO3--],[H+]) = 0 (=GG)
C ===========================================================
C
C
C* 2.1  SET MAX. NUMBER OF ITERATIONS
C --------------------------------------
C
      kitter = 30
C
C* 2.2  SET DAMPING PARAMETERS FOR CORRECTIONS OF [CO3--]
C       AND [H+]
C -------------------------------------------------------
C
C
C 2.3 INITIALISATION OF [HI+], and [CO3--]
C ----------------------------------------
C
      DO jk=1,jpkm1
        DO jj=1,jpj
          DO ji=1,jpi
            caralk = trn(ji,jj,jk,jptal)-
     &          borat(ji,jj,jk)/(1.+1E-8/akb3(ji,jj,jk))
            co3(ji,jj,jk) = caralk-trn(ji,jj,jk,jpdic)
     &          +(1.-tmask(ji,jj,jk))*.5e-3
            bicarb = (2.*trn(ji,jj,jk,jpdic)-caralk)
            hi(ji,jj,jk) = ak23(ji,jj,jk)*bicarb/co3(ji,jj,jk)
          END DO 
        END DO 
      END DO 
C
C* 2.4  BEGIN OF ITERATION
C ------------------------
C
      DO jn = 1,kitter
C
C* 2.5  COMPUTE [CO3--] and [H+] CONCENTRATIONS
C -------------------------------------------
C
        rconvs=0.
        DO jk = 1,jpkm1
          DO jj=1,jpj
            DO ji = 1, jpi
C
C* 2.6  SET DUMMY VARIABLE FOR TOTAL BORATE
C -----------------------------------------
C
              bot = borat(ji,jj,jk)
C
C* 2.7  SET DUMMY VARIABLE FOR [H+], AND [CO3--]
C ----------------------------------------------
C
              h = hi(ji,jj,jk)+(1.-tmask(ji,jj,jk))*1.e-9
              h = amax1(hi(ji,jj,jk),1.E-10)
              r = co3(ji,jj,jk)+(1.-tmask(ji,jj,jk))*.5e-3
C
C* 2.8  SET DUMMY VARIABLE FOR [ALK]GIVEN AND
C       [SUM(CO2)]GIVEN
C -------------------------------------------
C
              alka = trn(ji,jj,jk,jptal) 
              c = trn(ji,jj,jk,jpdic) 
C
C* 2.9 CALCULATE [ALK]([CO3--], [HCO3-])
C ------------------------------------
C
              a=alka-
     &            (akw3(ji,jj,jk)/h-h+bot/(1.+h/akb3(ji,jj,jk)))
C
C* 2.10 CALCULATE [H+] and [CO3--]
C -----------------------------------------
C
              ah2=sqrt((c-a)**2+4.*(a*ak23(ji,jj,jk)/ak13(ji,jj,jk))
     &            *(2*c-a))
              ah2=0.5*ak13(ji,jj,jk)/a*((c-a)+ah2)
              co3(ji,jj,jk) = a/(2.+ah2/ak23(ji,jj,jk))
C
C* 2.11 CONTROL VARIABLE TO CHECK CONVERGENCE
C -------------------------------------------
C
c
              conveg=((ah2-hi(ji,jj,jk))/hi(ji,jj,jk))**2
     $            *tmask(ji,jj,jk)
              rconvs = rconvs+conveg
              hi(ji,jj,jk)  = ah2
            ENDDO
          ENDDO
        END DO
C
C
C  2.12 CHECK CONVERGENCE
C  ----------------------
C
        IF (rconvs.LE.1.E-6) EXIT
C
      END DO 
C
C
C  2.13 CHECK ITERATION
C  --------------------
C
C      IF (jn.GE.kitter) THEN 
C          WRITE(numout,*) 'in h3clys hi and co3 did not converge '
C     $        ,rconvs
C          CALL FLUSH_(numout)
C          STOP 'h3clys'
C      END IF       
C
C     ---------------------------------------------------------
C*    3. CALCULATE DEGREE OF CACO3 SATURATION AND CORRESPONDING
C        DISSOLOUTION AND PRECIPITATION OF CACO3 (BE AWARE OF
C        MGCO3)
C     ---------------------------------------------------------
C
      DO jk = 2,jpkm1
        DO jj = 1,jpj
          DO ji = 1, jpi
C
C* 3.1  DEVIATION OF [CO3--] FROM SATURATION VALUE
C ------------------------------------------------
C
            delco3 = co3(ji,jj,jk)-aksp(ji,jj,jk)/calcon
C
C* 3.2  SET DEGREE OF UNDER-/SUPERSATURATION
C ------------------------------------------
C
            excess = amax1(zero,delco3)
C
C* 3.3  AMOUNT CACO3 (12C) THAT RE-ENTERS SOLUTION
C       (ACCORDING TO THIS FORMULATION ALSO SOME PARTICULATE
C       CACO3 GETS DISSOLVED EVEN IN THE CASE OF OVERSATURATION)
C --------------------------------------------------------------
C
            dispot = trn(ji,jj,jk,jpcal)*amin1(1.,
     &          (1.-delco3/(dispo0+abs(delco3))) )
#    if defined key_off_degrad
     &          *facvol(ji,jj,jk)
#    endif
C
C* 3.5  CHANGE OF PARTICULATE CACO3 AND TOTAL INORGANIC 14C
C       IN THE WATER COLUMN DUE TO CACO3 DISSOLUTION/PRECIP.
C ----------------------------------------------------------
C
            cristl = spocri
#    if defined key_off_degrad
     &          *facvol(ji,jj,jk)
#    endif
C* 3.8  CHANGE OF [CO3--] , [ALK], PARTICULATE [CACO3],
C       AND [SUM(CO2)] DUE TO CACO3 DISSOLUTION/PRECIPITATION
C -----------------------------------------------------------
C
            remco3=(dispot-excess*cristl)/rmoss
            co3(ji,jj,jk) = co3(ji,jj,jk)
     &          +remco3*rfact
            tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal)+
     &          2.*remco3
            tra(ji,jj,jk,jpcal) = tra(ji,jj,jk,jpcal)-
     &          remco3
            tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic)+
     &          remco3
#    if defined key_trc_biohamocc13
            tra(ji,jj,jk,jp13c) = tra(ji,jj,jk,jp13c)+pdb*
     &          remco3
#    endif
C
C
          ENDDO
        ENDDO
      END DO
C
#endif
      RETURN
      END

