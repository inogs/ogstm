      SUBROUTINE h3cexp
#if defined key_passivetrc && defined key_trc_hamocc3
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE h3cexp
CCC                     ******************
CCC
CC
CC     PURPOSE.
CC     --------
CC          *H3CEXP* MODELS EXPORT OF BIOGENIC MATTER (POC ''SOFT
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
CC          INTERNATIONAL ASSOCIATION OF METEOROLOGY AND ATMOSPHERIC
c PHYSICS.
CC
CC   MODIFICATIONS:
CC   --------------
CC      original      : 1999    O. Aumont
CC      modifications : 1999    C. Le Quere
CC                    : 2001    O. Aumont
CC                    : 01/03   O. Aumont,EK: add sediments 
CC                              for calcite and silicate
CC ---------------------------------------------------------------------
c ------
CC parameters and commons
CC ======================


       USE myalloc
       USE stafun

CC----------------------------------------------------------------------
CC local declarations
CC ==================
C
      INTEGER ji, jj, jk, zkbot(jpi,jpj)
      REAL(8) prcaca12(jpi,jpj),silpro12(jpi,jpj)
      REAL(8) zarea, zvol,silcri,sister
      REAL(8) zgeolsil, zgeolcal
      REAL(8) sumcal
C
      silcri = 4.e-4*rfact/(3600.*24.*30.5)
      prcaca12 = 0.0
      fbod = 0.0
      silpro12 = 0.0
C
C VERTICAL DISTRIBUTION OF NEWLY PRODUCED BIOGENIC
C CACO3 PARTICLES AND POC IN THE WATER COLUMN
C (PARTS OF NEWLY FORMED MATTER REMAINING IN THE DIFFERENT
C LAYERS IS DETERMINED BY DMIN3, DISC3 AS DEFINED IN common.passivetrc
c .hamocc3.h
C ----------------------------------------------------------------------
c --------
C
      DO jk = 1,jpkb
        DO jj = 1,jpj
          DO ji = 1,jpi
            prcaca12(ji,jj) = prcaca12(ji,jj)+
     &          e3t(jk)*prcaca(ji,jj,jk)*tmask(ji,jj,jk)*rfactr
            silpro12(ji,jj) = silpro12(ji,jj)+silpro(ji,jj,jk)*
     &          tmask(ji,jj,jk)*e3t(jk)*rfactr
          ENDDO
        ENDDO
      ENDDO

C
#    if defined key_trc_p3zd
      DO jj = 1,jpj
        DO ji = 1,jpi
          fbod(ji,jj) = fbod(ji,jj)+(sinking(ji,jj,11)
     $        +nu(ji,jj,11))*tmask(ji,jj,10)*rfactr
        ENDDO
      ENDDO
#    else
      DO jk = 1,jpkb
        DO jj = 1,jpj
          DO ji = 1,jpi
            fbod(ji,jj) = fbod(ji,jj)+
     &          e3t(jk)*prorca(ji,jj,jk)*tmask(ji,jj,jk)
     &          *rfactr
          ENDDO
        ENDDO
      ENDDO
#    endif
C
#    if defined key_trc_diaadd
      DO jj = 1,jpj
        DO ji = 1,jpi
          trc2d(ji,jj,5) = fbod(ji,jj)*1.e3
          trc2d(ji,jj,6) = prcaca12(ji,jj)*1.e3
          trc2d(ji,jj,7) = silpro12(ji,jj)*1.e3
        ENDDO
      ENDDO 
#    endif
C
      DO jk = 1,jpk
        DO jj = 1,jpj
          DO ji = 1,jpi
#    if defined key_trc_p3zd
            trn(ji,jj,jk,jpcal) = trn(ji,jj,jk,jpcal)+disc3(ji,jj,jk)*
     &          prcaca12(ji,jj)/e3t(jk)*rdt
            trn(ji,jj,jk,jpsil) = trn(ji,jj,jk,jpsil)+silpro12(ji,jj)*
     &          diss3(ji,jj,jk)*1./e3t(jk)*rdt
CC          trn(ji,jj,jk,jpsil) = min(200.e-6,trn(ji,jj,jk,jpsil))
#    else
            tra(ji,jj,jk,jpcal) = tra(ji,jj,jk,jpcal)+disc3(ji,jj,jk)*
     &                            prcaca12(ji,jj)/e3t(jk)
            tra(ji,jj,jk,jpsil) = tra(ji,jj,jk,jpsil)+silpro12(ji,jj)*
     &                            diss3(ji,jj,jk)*1./e3t(jk)
CC          sister = max(0.,silpro12(ji,jj)*diss3(ji,jj,jk)-silcri)
CC          tra(ji,jj,jk,jpsil) = tra(ji,jj,jk,jpsil)-
CC   &                            sister/e3t(jk)*rfactr
#    endif
          ENDDO
        ENDDO
      ENDDO
C
      zkbot = jpk
C
      DO jk = 1,jpkm1
        DO jj = 2,jpjm1
          DO ji = 2,jpim1
  
             IF (tmask(ji,jj,jk).eq.1.and.
     .           tmask(ji,jj,jk+1).eq.0) THEN
C
                  zkbot(ji,jj) = jk
C
              ENDIF

             ENDDO
         ENDDO
      ENDDO
C
C
C     Initialization
      zgeolcal = 0.
      zgeolsil = 0.

C
C     sedlam --> sedlam/100 for silicate and calcite
C
        DO jj = 2,jpjm1
          DO ji = 2,jpim1
#            if defined key_trc_p3zd
             trn(ji,jj,zkbot(ji,jj),jpcal) =
     .          trn(ji,jj,zkbot(ji,jj),jpcal) +
     .               sedlam/100.
     .               *sedcal(ji,jj)*rdt/fse3t(ji,jj,zkbot(ji,jj))
             trn(ji,jj,zkbot(ji,jj),jpsil) =
     .          trn(ji,jj,zkbot(ji,jj),jpsil) +
     .               sedlam/100.
     .               *sedsil(ji,jj)*rdt/fse3t(ji,jj,zkbot(ji,jj))
#            elif defined key_trc_hamocc3 && ! defined key_trc_p3zd
             tra(ji,jj,zkbot(ji,jj),jpcal) =
     .          tra(ji,jj,zkbot(ji,jj),jpcal) +
     .               sedlam/100.
     .               *sedcal(ji,jj)/fse3t(ji,jj,zkbot(ji,jj))
             tra(ji,jj,zkbot(ji,jj),jpsil) =
     .          tra(ji,jj,zkbot(ji,jj),jpsil) +
     .               sedlam/100.
     .               *sedsil(ji,jj)/fse3t(ji,jj,zkbot(ji,jj))
#            endif

             zgeolcal = zgeolcal + sedlostcal*sedcal(ji,jj)*
     .                             e1t(ji,jj)*e2t(ji,jj)
             sedcal(ji,jj) = sedcal(ji,jj)  +
     .                       discl(ji,jj)*prcaca12(ji,jj)*rdt -
     .                       sedlam/100.*sedcal(ji,jj)*rdt -
     .                       sedlostcal*sedcal(ji,jj)*rdt

             zgeolsil = zgeolsil + sedlostsil*sedsil(ji,jj)*
     .                             e1t(ji,jj)*e2t(ji,jj)
             sedsil(ji,jj) = sedsil(ji,jj)  +
     .                       dissl(ji,jj)*silpro12(ji,jj)*rdt -
     .                       sedlam/100.*sedsil(ji,jj)*rdt - 
     .                       sedlostsil*sedsil(ji,jj)*rdt

C
             ENDDO
         ENDDO

C     sumcal = 0.
C      DO jj = 1,jpj
C         DO ji = 1,jpi
C         sumcal = sumcal + tmask(ji,jj,1)* 
C    .         sedcal(ji,jj) * e1t(ji,jj)*e2t(ji,jj)
C       ENDDO
C     ENDDO

C     write(0,*) 'sumcal2 ',sumcal*1000.*12/1E15
C
C
CC 12/1/2004 Paolo - casi 0/0 ================================
      if (areacot.ne.0) then

        DO jj = 2,jpjm1
          DO ji = 2,jpim1

#            if defined key_trc_p3zd
             trn(ji,jj,1,jpcal) = trn(ji,jj,1,jpcal) + zgeolcal*rdt*
     .                            cmask(ji,jj)/areacot/fse3t(ji,jj,1)
             trn(ji,jj,1,jpsil) = trn(ji,jj,1,jpsil) + zgeolsil*rdt*
     .                            cmask(ji,jj)/areacot/fse3t(ji,jj,1)
#            elif defined key_trc_hamocc3 && ! defined key_trc_p3zd
             tra(ji,jj,1,jpcal) = tra(ji,jj,1,jpcal) + zgeolcal*
     .                            cmask(ji,jj)/areacot/fse3t(ji,jj,1)
             tra(ji,jj,1,jpsil) = tra(ji,jj,1,jpsil) + zgeolsil*
     .                            cmask(ji,jj)/areacot/fse3t(ji,jj,1)
#            endif

             ENDDO
         ENDDO
      endif
CC ==========================================================
C
         CALL lbc(sedcal,1,1,1,1,1,1)
         CALL lbc(sedsil,1,1,1,1,1,1)

C Oa & Ek: diagnostics depending on jpdia2d
C          left as example
#     if defined key_trc_diaadd
           do jj=1,jpj
             do ji=1,jpi
              trc2d(ji,jj,12)=sedcal(ji,jj)
              trc2d(ji,jj,13)=sedsil(ji,jj)
             end do
           end do
#     endif

#endif

      RETURN
      END
