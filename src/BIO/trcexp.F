C $Id: 
      SUBROUTINE trcexp
#if defined key_passivetrc 
#if defined key_trc_npzd || defined key_trc_nnpzddom || defined key_trc_hamocc3
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE trcexp
CCC                     ******************
CCC
CC
CC     PURPOSE.
CC     --------
CC          *TRCEXP* MODELS EXPORT OF BIOGENIC MATTER (POC ''SOFT
CC                   TISSUE'') AND ITS DISTRIBUTION IN WATER COLUMN
CC
CC     METHOD.
CC     -------
CC          IN THE SURFACE LAYER POC IS PRODUCED ACCORDING TO
CC     NURTRIENTS AVAILABLE AND GROWTH CONDITIONS. NUTRIENT UPTAKE
CC     KINETICS FOLLOW MICHAELIS-MENTON FORMULATION. 
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
CC   MODIFICATIONS:
CC   --------------
CC      original      : 1999    O. Aumont
CC      modifications : 1999    C. Le Quere
CC      additions   : 01-05 (O. Aumont, E. Kestenare):
CC                           add sediment computations
CC ---------------------------------------------------------------------

       USE myalloc
       USE stafun
        IMPLICIT NONE


CC----------------------------------------------------------------------
CC local declarations
CC ==================
C
      INTEGER ji, jj, jk, zkbot(jpi,jpj)
      REAL(8) zwork(jpi,jpj), zgeolpoc


C
C VERTICAL DISTRIBUTION OF NEWLY PRODUCED BIOGENIC
C POC IN THE WATER COLUMN
C (PARTS OF NEWLY FORMED MATTER REMAINING IN THE DIFFERENT
C LAYERS IS DETERMINED BY DMIN3 DEFINED IN common.passivetrc.*.h
C ----------------------------------------------------------------------
C
C
      DO jk = 1,jpkm1
        DO jj = 2,jpjm1
          DO ji = 2,jpim1
#    if defined key_trc_p3zd
            trn(ji,jj,jk,jppoc) = trn(ji,jj,jk,jppoc)+
     &          (1./fse3t(ji,jj,jk))*rdt*
     &          dmin3(ji,jj,jk) *fbod(ji,jj)
#    elif defined key_trc_hamocc3 && ! defined key_trc_p3zd
            tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc)+
CCC 28/4/2004 Paolo time step term rdt*
CCC     &          (1./fse3t(ji,jj,jk))*
     &          (1./fse3t(ji,jj,jk))*rdt*
     &          dmin3(ji,jj,jk) *fbod(ji,jj)
#    else
            tra(ji,jj,jk,jpno3) = tra(ji,jj,jk,jpno3)+
CCC 28/4/2004 Paolo time step term rdt*
CCC     &          (1./fse3t(ji,jj,jk))*
     &          (1./fse3t(ji,jj,jk))*rdt*
     &          dmin3(ji,jj,jk) *fbod(ji,jj)
#    endif
          ENDDO
        ENDDO
      ENDDO
C
C     Find the last level of the water column
C     Compute fluxes due to sinking particles (slow)
C   
      zkbot = jpk
      zwork = 0.
C
C
      DO jk = 1,jpkm1
        DO jj = 2,jpjm1
          DO ji = 2,jpim1
   
             IF (tmask(ji,jj,jk).eq.1.and.
     .           tmask(ji,jj,jk+1).eq.0) THEN
C
                  zkbot(ji,jj) = jk
#    if  ! defined key_trc_hamocc3
CCC 26/4/2004 Paolo speed from (m/day) to (m/s) 
                  zwork(ji,jj) = vsed * trn(ji,jj,jk,jpdet)/24./3600.
#    endif
C
              ENDIF
    
             ENDDO
         ENDDO
      ENDDO
C
C     Initialization
      zgeolpoc = 0.

C     Release of nutrients from the "simple" sediment
C
        DO jj = 2,jpjm1
          DO ji = 2,jpim1

#            if defined key_trc_p3zd
             trn(ji,jj,zkbot(ji,jj),jppo4) = 
     .          trn(ji,jj,zkbot(ji,jj),jppo4) +
     .               sedlam*sedpoc(ji,jj)*rdt/fse3t(ji,jj,zkbot(ji,jj))
#            elif defined key_trc_hamocc3 && ! defined key_trc_p3zd
             tra(ji,jj,zkbot(ji,jj),jppo4) = 
     .          tra(ji,jj,zkbot(ji,jj),jppo4) +
     .               sedlam*sedpoc(ji,jj)/fse3t(ji,jj,zkbot(ji,jj))
#            else
CCC 26/4/2004 Paolo time step(*rdt) 
             tra(ji,jj,zkbot(ji,jj),jpno3) = 
     .          tra(ji,jj,zkbot(ji,jj),jpno3) +
     .               sedlam*sedpoc(ji,jj)/fse3t(ji,jj,zkbot(ji,jj))*rdt
#            endif
C
C     Deposition of organic matter in the sediment
C
             zgeolpoc = zgeolpoc + sedlostpoc*sedpoc(ji,jj)*
     .                             e1t(ji,jj)*e2t(ji,jj)

             sedpoc(ji,jj) = sedpoc(ji,jj)  +
     .                       zwork(ji,jj)*rdt +
     .                       dminl(ji,jj)*fbod(ji,jj)*rdt -
     .                       sedlam*sedpoc(ji,jj)*rdt -
     .                       sedlostpoc*sedpoc(ji,jj)*rdt
C
             ENDDO
         ENDDO
CC 12/1/2004 Paolo - casi 0/0 ================================
      if (areacot.ne.0) then

C
        DO jj = 2,jpjm1
          DO ji = 2,jpim1

#            if defined key_trc_p3zd
             trn(ji,jj,1,jppo4) = trn(ji,jj,1,jppo4) + zgeolpoc*rdt*
     .                            cmask(ji,jj)/areacot/fse3t(ji,jj,1)
#            elif defined key_trc_hamocc3 && ! defined key_trc_p3zd
             tra(ji,jj,1,jppo4) = tra(ji,jj,1,jppo4) + zgeolpoc*
     .                            cmask(ji,jj)/areacot/fse3t(ji,jj,1)
#            else
CCC 26/4/2004 Paolo time step(*rdt)
             tra(ji,jj,1,jpno3) = tra(ji,jj,1,jpno3) + zgeolpoc*
     .                            cmask(ji,jj)/areacot/fse3t(ji,jj,1)*rdt
#            endif

             ENDDO
         ENDDO
      endif
CC ==========================================================

         CALL lbc(sedpoc,1,1,1,1,1,1)

C Oa & Ek: diagnostics depending on jpdia2d
C          left as example
#     if defined key_trc_diaadd
           do jj=1,jpj
             do ji=1,jpi
              trc2d(ji,jj,11)=sedpoc(ji,jj)
C              trc2d(ji,jj,5) = fbod(ji,jj)
             end do
           end do
#     endif

#        if defined key_trc_p3zd
            CALL lbc(trn,1,1,1,1,1,1)
#        endif
C
#endif
#endif

      END SUBROUTINE trcexp
