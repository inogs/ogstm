      SUBROUTINE h3crem
#if defined key_passivetrc && defined key_trc_hamocc3
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE h3crem
CCC                     ******************
CCC
CCC
CCC
CCC
CCC     PURPOSE.
CCC     --------
CCC          *H3CREM* CALCULATES OXYGEN CONSUMPTION DURING REMINERALIZATION
CCC                   AND LOSS OF POC TO THE ORGANIC SEDIMENT POOL
CCC
CC
CC     METHOD.
CC     -------
CC          THE OXYGEN AMOUNT AVAILABLE FOR CONSUMPTION DUE TO
CC     REMINERALIZATION OF POC IS DETERMINED.
CC          THE CORRESPONDING DECREASE OF POC AND THE INCREASES OF
CC     INORGANIC CARBON AND PHOSPHATE ARE CALCULATED. ALKALINITY
CC     IS REDUCED BY THE NITRATE SET FREE DURING REMINERALIZATION.
CC          POC SETTLES DOWN TO THE SEDIMENT FROM THE WATER LAYER
CC     DIRECTLY ABOVE THE OCEAN FLOOR (DEPENDS ON BOTTOM DEPTH)
CC     AT A RATE TO BE DEFINED (''OCSED'').
CC
CC     EXTERNALS.
CC     ----------
CC          NONE.
CC
CC     REFERENCE.
CC     ----------
CC          -.
CC
CC*     VARIABLE      TYPE    PURPOSE.
CC      --------      ----    --------
CC
CC      *ZREM*        REAL(8)    MOLES O2 AVAILABLE FOR REMINERALIZATION
CC                            (LIMITING INGREDIENT FOR REMINERALIZATION)
CC      *ZREMIK*      REAL(8)    REMINERALIZATION RATE (FRACTION OF POC,
CC                            THAT IS REMINERALISED DURING ONE TIMESTEP)
CC                            DUMMY VARIABLE
CC
CC   MODIFICATIONS:
CC   --------------
CC      original      : 1988-07 E. MAIER-REIMER      MPI HAMBURG
CC      additions     : 1998    O. Aumont
CC      modifications : 1999    C. Le Quere
CC ---------------------------------------------------------------------
CC  parameters and commons
CC ======================

       USE myalloc
       IMPLICIT NONE

CC----------------------------------------------------------------------
CC local declarations
CC ==================
C
      INTEGER ji, jj, jk
      REAL(8) zremik, zrem
C
C* 3. REMINERALIZATION OF POC AND O2 UTILIZATION
C ----------------------------------------------
C
C* 3.1  LAYER LOOP
C ----------------
C
      DO jk = 1,jpkm1
        DO jj = 1,jpj
C
C* 3.2  SET FRACTION POC TO BE REMINERALISED IN LAYER K
C -----------------------------------------------------
C
C* 3.3  MERIDIONAL/ZONAL LOOPS
C ----------------------------
C
          DO ji = 1,jpi
C
            zrem = 0.
C
            zremik = 50.*spocri*(1.-tmask(ji,jj,jk+1))
     $          +0.002*(tn(ji,jj,jk)+2.)
            zremik = 0.002*(tn(ji,jj,jk)+2.)
            zremik = zremik/(30.5*3600.*24.)*rfact
# if defined key_off_degrad
     &          *facvol(ji,jj,jk)
# endif
C
C* 3.5  AMOUNT OF O2 AVAILABLE FOR REMINERALIZATION OF POC
C --------------------------------------------------------
C
            zrem = zremik* min(
     &          (trn(ji,jj,jk,jpoxy)-oxymin)/o2ut,trn(ji,jj,jk,jppoc))
C
C
C* 3.6  INCREASE OF PO4 AND SUM((12C)O2) DUE TO RELEASE TO WATER
C       COLUMN DURING REMINERALIZATION
C --------------------------------------------------------------
C
            tra(ji,jj,jk,jppo4) = tra(ji,jj,jk,jppo4)+ zrem*rfactr
            tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic)+ zrem*rfactr
C
C* 3.7  DECREASE OF ALKALINITY AND O2 DUE TO REMINERALIZATION
C       (ALKALINITY DECREASES DUE TO NITRATE RELEASE)
C -----------------------------------------------------------
C
            tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal)-
     &          zrem*rno3*rfactr
            tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy)- 
     &          rfactr*o2ut*zrem
C
C* 3.8  CHANGE IN RARE ISOTOPES OF INORGANIC CARBON AND POC
C       AFTER REMINERALIZATION AND DECREASE OF POC 12C
C       IN THE WATER COLUMN DUE TO REMINERALIZATION
C ----------------------------------------------------------------------
C
#    if defined key_trc_biohamocc13
            tra(ji,jj,jk,jp13c) = tra(ji,jj,jk,jp13c)+
     &          zrem*pdb*rfactr*plafr13
#    endif
            tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc)-
     &          zrem*rfactr
C
          ENDDO
        ENDDO
C
C END of slab
C ===========
C
      ENDDO 
# endif
      RETURN
      END
      
