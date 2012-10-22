C $Id: p3zrem.F,v 1.2 2009-09-11 09:20:56 cvsogs01 Exp $
       SUBROUTINE p3zrem
#if defined key_passivetrc && key_trc_p3zd
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE p3zrem
CCC                     ******************
CCC
CCC
CCC
CCC
CCC     PURPOSE.
CCC     --------
CC
CC         *P3ZREM* CALCULATES OXYGEN CONSUMPTION DURING REMINERALIZATION
CC                  AND LOSS OF POC TO THE ORGANIC SEDIMENT POOL
CC
CC*   INTERFACE.
CC    ----------
CC         *CALL* *P3ZREM*
CC
CC    METHOD.
CC    -------
CC         THE OXYGEN AMOUNT AVAILABLE FOR CONSUMPTION DUE TO
CC    REMINERALIZATION OF POC IS DETERMINED.
CC         THE CORRESPONDING DECREASE OF POC AND THE INCREASES OF
CC    INORGANIC CARBON AND PHOSPHATE ARE CALCULATED. ALKALINITY
CC    IS REDUCED BY THE NITRATE SET FREE DURING REMINERALIZATION.
CC         POC SETTLES DOWN TO THE SEDIMENT FROM THE WATER LAYER
CC    DIRECTLY ABOVE THE OCEAN FLOOR (DEPENDS ON BOTTOM DEPTH)
CC    AT A RATE TO BE DEFINED (''OCSED'').
CC
CC    EXTERNALS.
CC    ----------
CC         NONE.
CC
CC    REFERENCE.
CC    ----------
CC         -.
CC
CC     VARIABLE      TYPE    PURPOSE.
CC     --------      ----    --------
CC
CC     *OLIMI(IE,JE)*REAL(8)    MOLES O2 AVAILABLE FOR REMINERALIZATION
CC                           (LIMITING INGREDIENT FOR REMINERALIZATION)
CC     *KE1*         INTEGER MAX. OF INDEX FOR OCEAN LAYERS MINUS 1
CC     *OCSED*       REAL(8)    FRACTION OF POC LOST TO ORGANIC CARBON SEDIMENT
CC                           POOL FROM WATER LAYER IMMEDIATELY ABOVE OCEAN
CC                           FLOOR
CC     *OCSEDI*      REAL(8)    FRACTION OF POC REMAINING IN THE WATER LAYER
CC                           IMMEDIATELY ABOVE THE BOTTOM AND DOES NOT
CC                           ENTER THE ORG. CARBON SEDIMENT POOL
CC     *REMIK*       REAL(8)    REMINERALIZATION RATE (FRACTION OF POC,
CC                           THAT IS REMINERALISED DURING ONE TIMESTEP),
CC                           DUMMY VARIABLE
CC
CC   MODIFICATIONS:
CC   --------------
CC      original      : 1988-07 E. MAIER-REIMER      MPI HAMBURG
CC      additions     : 1998    O. Aumont
CC      modifications : 1999    C. Le Quere
CC ----------------------------------------------------------------------
CC  parameters and commons
CC ======================

       USE myalloc

        IMPLICIT NONE


CC----------------------------------------------------------------------
CC local declarations
CC ==================
C
      INTEGER ji, jj, jk, jn, ke1
      INTEGER ktask
      REAL(8) zrem, zorem, zremip, zremik, zphy, zzoo
C
C
C REMINERALIZATION OF POC AND O2 UTILIZATION
C ------------------------------------------
C
      DO jk = 11,jpk-1
C
C SET FRACTION POC TO BE REMINERALISED IN LAYER K
C -----------------------------------------------
C
        DO jj = 1,jpj
          DO ji = 1,jpi
C
            zremik=50.*spocri*(1.-tmask(ji,jj,jk+1))+xremik
            zremik = zremik/rjjss*rfact*tmask(ji,jj,jk)
     &          *trn(ji,jj,jk,jppo4)/(trn(ji,jj,jk,jppo4)+xkdoc1)
     &          *trn(ji,jj,jk,jpdoc)/(trn(ji,jj,jk,jpdoc)+xkdoc2)
#    if defined key_off_degrad
     &          *facvol(ji,jj,jk)
#    endif
            zremip = xremip/rjjss*rfact*tmask(ji,jj,jk)
#    if defined key_off_degrad
     &          *facvol(ji,jj,jk)
#    endif
            zzoo = (resrat+mzrat)/rjjss*rfact*tmask(ji,jj,jk)*
     &          trn(ji,jj,jk,jpzoo)*trn(ji,jj,jk,jpzoo)/
     &          (1.E-6+trn(ji,jj,jk,jpzoo))
#    if defined key_off_degrad
     &          *facvol(ji,jj,jk)
#    endif
            zphy = mprat/rjjss*rfact*tmask(ji,jj,jk)*
     &          trn(ji,jj,jk,jpphy)
     &          *trn(ji,jj,jk,jpphy)/(1.E-6+trn(ji,jj,jk,jpphy))
#    if defined key_off_degrad
     &          *facvol(ji,jj,jk)
#    endif
C
C AMOUNT OF O2 AVAILABLE FOR REMINERALIZATION OF POC
C
            zrem = zremik* trn(ji,jj,jk,jpdoc) 
            zorem = zremip*trn(ji,jj,jk,jppoc)
C
C
C INCREASE OF PO4 AND SUM((12C)O2) DUE TO RELEASE TO WATER
C COLUMN DURING REMINERALIZATION
C --------------------------------------------------------
C
            trn(ji,jj,jk,jpzoo) = trn(ji,jj,jk,jpzoo)-zzoo
            trn(ji,jj,jk,jpphy) = trn(ji,jj,jk,jpphy)-zphy
            trn(ji,jj,jk,jppo4) = trn(ji,jj,jk,jppo4)+zrem
            trn(ji,jj,jk,jpdic) = trn(ji,jj,jk,jpdic)+zrem
C
C DECREASE OF ALKALINITY AND O2 DUE TO REMINERALIZATION
C (ALKALINITY DECREASES DUE TO NITRATE RELEASE)
C -----------------------------------------------------
C
            trn(ji,jj,jk,jptal) = trn(ji,jj,jk,jptal)-zrem*rno3
            trn(ji,jj,jk,jpoxy) = trn(ji,jj,jk,jpoxy)-o2ut*zrem
            trn(ji,jj,jk,jpoxy) = MAX(4.E-6,trn(ji,jj,jk,jpoxy))
C
#    if defined key_trc_biohamocc13
C CHANGE IN RARE ISOTOPES OF INORGANIC CARBON AND POC
C AFTER REMINERALIZATION AND DECREASE OF POC 12C
C IN THE WATER COLUMN DUE TO REMINERALIZATION
C
            trn(ji,jj,jk,jp13c) = trn(ji,jj,jk,jp13c)+zrem*pdb*
     &                           plafr13
#    endif
C
            trn(ji,jj,jk,jppoc) = trn(ji,jj,jk,jppoc)-zorem
            trn(ji,jj,jk,jpdoc) = trn(ji,jj,jk,jpdoc)+zorem+zphy
     &                           +zzoo-zrem
C
          ENDDO
        ENDDO
      ENDDO
C
#endif
C
      RETURN
      END
