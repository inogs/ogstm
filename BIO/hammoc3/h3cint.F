C $Id: h3cint.F,v 1.2 2009-09-11 09:20:56 cvsogs01 Exp $

      SUBROUTINE h3cint(kt)
C
#if defined key_trc_hamocc3
C

c#include "parameter.h"
c#include "common.h"
       USE myalloc

CC+CC  Implicit typing is never allowed
        IMPLICIT NONE
CC+CC  Implicit typing is never allowed

      INTEGER kt
      INTEGER ji, jj, jk, l, imel(jpi,jpj)
      INTEGER iman
      INTEGER iv, ipdtant, ipdtmo
      INTEGER nvit1t,nvit2t
      REAL(8)       zpdtan, zman, zpdtmo, zdemi
      REAL(8) xtxy, zt
      REAL(8) zavmx
C
      ipdtant  = (12*730*60*60)/int(rdttra(1))
      iman    = ilecvit
      ipdtmo  = ipdtant / iman
C
      zpdtan = float ( ipdtant)
      zman   = float ( iman)
      zpdtmo = zpdtan / zman
      zdemi  = zpdtmo / 2.
      zt     = ( float ( kt) + zdemi) / zpdtmo
      
C  recherche de l'indice des enregistrements
C  du modele dynamique encadrant le pas de temps kt.
C  --------------------------------------------------
C
      xtvit = zt - int ( zt)
      nvit1t = int (( float ( kt) + zdemi)/ zpdtmo)
      nvit2t = nvit1t+1
      nvit1t = MOD ( nvit1t, iman)
      IF ( nvit1t .EQ. 0 ) nvit1t = iman
      nvit2t = MOD ( nvit2t, iman)
      IF ( nvit2t .EQ. 0 ) nvit2t = iman
C
C vertical slab
C ===============
C
      IF (icice.NE.3) THEN
          DO jj = 1,jpj
C
c  ---------------------------------------------------------
c  interpolation de la couverture de glace et du vent
c  ---------------------------------------------------------
c
            DO ji=1,jpi
              freeze(ji,jj) = (1.-xtvit)*cicemo(ji,jj,nvit1t)
     $            +xtvit*cicemo(ji,jj,nvit2t)
            END DO
          END DO
      ENDIF
      IF (igaswind.NE.3) THEN
          DO jj = 1,jpj
            DO ji=1,jpi
              vatm(ji,jj) = (1.-xtvit)*wsmo(ji,jj,nvit1t)
     $            +xtvit*wsmo(ji,jj,nvit2t)
            ENDDO
          END DO
      ENDIF
      IF (igaswind.EQ.2) THEN
          DO jj = 1,jpj
            DO ji=1,jpi
              kgwanin(ji,jj) = (1.-xtvit)*kgwanmo(ji,jj,nvit1t)
     $            +xtvit*kgwanmo(ji,jj,nvit2t)
            ENDDO
          END DO
      ENDIF
C
C
#endif
C
      RETURN
      END


