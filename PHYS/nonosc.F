
      SUBROUTINE nonosc(pbef,paa,pbb,pcc,paft,prdt)
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE nonosc
CCC                      ****************
CCC
CCC  PURPOSE :
CCC  ---------
CCC     compute monotonic tracer fluxes from the upstream scheme
CCC     and the before field by a nonoscillatory algorithm 
CCC
CC   METHOD :
CC   -------
CC      warning : pbef and paft must be masked, but the boundaries
CC                conditions on the fluxes are not necessary
CC      m.zalezak (1979)
CC      h.drange (1995) multi-dimensional forward-in-time and
CC                      upstream-in-space based differencing
CC                      for fluid
CC
CC   INPUT :
CC   -----
CC      argument
CC              pbef            : before field (tracer)
CC              paa             : mass flux per second (direction i)
CC              pbb             : mass flux per second (direction j)
CC              pcc             : mass flux per second (direction k)
CC              paft            : after field (tracer) (upstrem scheme)
CC              prdt            : time step
CC
CC   OUTPUT :
CC   ------
CC              paa             : mass flux per second (direction i)
CC              pbb             : mass flux per second (direction j)
CC              pcc             : mass flux per second (direction k)
CC
CC   WORKSPACE :
CC   ---------
CC      local
CC              zkx(),zky(),zkz(),zbetup(),zbetdo()
CC      COMMON
CC
CC   EXTERNAL :   bord3
CC   --------
CC
CC   MODIFICATIONS:
CC   --------------
CC      original : 97-04 (l. mortier)
CC      addition : 00-02 (h. loukos) rewritting for opa8
CC      addition : 00-10 (M.A Foujols, E. Kestenare) 
CC                       mpplnk and lbc instead of bord3
CC      addition : 01-03 (E. Kestenare) add key_passivetrc
CC
CC----------------------------------------------------------------------
CC----------------------------------------------------------------------
CC parameters and commons
CC ======================

c#include "parameter.h"
c#include "common.h"
       USE myalloc

CC+CC  Implicit typing is never allowed
        IMPLICIT NONE
CC+CC  Implicit typing is never allowed

CC----------------------------------------------------------------------
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      REAL(8) prdt
      REAL(8) pbef(jpi,jpj,jpk),paft(jpi,jpj,jpk)
      REAL(8) paa(jpi,jpj,jpk),pbb(jpi,jpj,jpk),pcc(jpi,jpj,jpk)
#if defined key_passivetrc
      INTEGER ji,jj,jk
      INTEGER ikm1
      REAL(8) zbetup(jpi,jpj,jpk),zbetdo(jpi,jpj,jpk)
      REAL(8) zpos, zneg, zbt, za, zb, zc, zbig
CC----------------------------------------------------------------------
CC statement functions
CC ===================

#include "stafun.h"

CCC---------------------------------------------------------------------
CCC  opa8, LODYC (february 2000)
CCC---------------------------------------------------------------------
C
      zbig = 1.e+40
C
C
C 1. Search local extrema
C -----------------------
C to account for boundaries which are coded with tmask=0
      DO jk=1,jpk
        DO jj=1,jpj
          DO ji=1,jpi
            IF (tmask(ji,jj,jk).EQ.0.) THEN
                pbef(ji,jj,jk) = -zbig
                paft(ji,jj,jk) = -zbig
            ENDIF
          END DO
        END DO
      END DO
C search maximum in neighbourhood
      DO jk=1,jpkm1
        ikm1=MAX(jk-1,1)
        DO jj=2,jpjm1
          DO ji=2,jpim1
            zbetup(ji,jj,jk)=MAX(
     &          pbef(ji,jj,jk),  paft(ji,jj,jk),
     &          pbef(ji-1,jj,jk),pbef(ji+1,jj,jk),
     &          paft(ji-1,jj,jk),paft(ji+1,jj,jk),
     &          pbef(ji,jj-1,jk),pbef(ji,jj+1,jk),
     &          paft(ji,jj-1,jk),paft(ji,jj+1,jk),
     &          pbef(ji,jj,ikm1),pbef(ji,jj,jk+1),
     &          paft(ji,jj,ikm1),paft(ji,jj,jk+1))
          END DO
        END DO
      END DO
C to account for boundaries which are coded with tmask=0
      DO jk=1,jpk
        DO jj=1,jpj
          DO ji=1,jpi
            IF (tmask(ji,jj,jk).EQ.0.) THEN
                pbef(ji,jj,jk) = +zbig
                paft(ji,jj,jk) = +zbig
            ENDIF
          END DO
        END DO
      END DO
C search minimum in neighbourhood
      DO jk=1,jpkm1
        ikm1=MAX(jk-1,1)
        DO jj=2,jpjm1
          DO ji=2,jpim1
            zbetdo(ji,jj,jk)=MIN(
     &          pbef(ji,jj,jk),  paft(ji,jj,jk),
     &          pbef(ji-1,jj,jk),pbef(ji+1,jj,jk),
     &          paft(ji-1,jj,jk),paft(ji+1,jj,jk),
     &          pbef(ji,jj-1,jk),pbef(ji,jj+1,jk),
     &          paft(ji,jj-1,jk),paft(ji,jj+1,jk),
     &          pbef(ji,jj,ikm1),pbef(ji,jj,jk+1),
     &          paft(ji,jj,ikm1),paft(ji,jj,jk+1))
          END DO
        END DO
      END DO
C
      DO jk=1,jpk
        DO jj=1,jpj
          DO ji=1,jpi
C restore masked values
            IF (tmask(ji,jj,jk).EQ.0.) THEN
                pbef(ji,jj,jk) = 0.
                paft(ji,jj,jk) = 0.
            ENDIF
          END DO
        END DO
      END DO
C 
C
C 2. Positive and negative part of fluxes and beta terms
C ------------------------------------------------------
C
      DO jk=1,jpkm1
        DO jj=2,jpjm1
          DO ji=2,jpim1
            zpos= MAX(0.D0,paa(ji-1,jj,jk))-MIN(0.D0,paa(ji,jj,jk))
     &           +MAX(0.D0,pbb(ji,jj-1,jk))-MIN(0.D0,pbb(ji,jj,jk))
     &           +MAX(0.D0,pcc(ji,jj,jk+1))-MIN(0.D0,pcc(ji,jj,jk))
            zneg= MAX(0.D0,paa(ji,jj,jk))-MIN(0.D0,paa(ji-1,jj,jk))
     &           +MAX(0.D0,pbb(ji,jj,jk))-MIN(0.D0,pbb(ji,jj-1,jk))
     &           +MAX(0.D0,pcc(ji,jj,jk))-MIN(0.D0,pcc(ji,jj,jk+1))
C beta terms up and down
            zbt= e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk)/prdt
            zbetup(ji,jj,jk)=(zbetup(ji,jj,jk)-paft(ji,jj,jk))
     &          /(zpos+rtrn)*zbt
            zbetdo(ji,jj,jk)=(paft(ji,jj,jk)-zbetdo(ji,jj,jk))
     &          /(zneg+rtrn)*zbt
          END DO
        END DO
      END DO
C
Cek: 27/10/2000: call mpplnk or lbc instead of bord3
Cek   CALL bord3(zbetup)
Cek   CALL bord3(zbetdo)
#     ifdef key_mpp
        CALL mpplnk( zbetup(1,1,1), 1, 1 )
        CALL mpplnk( zbetdo(1,1,1), 1, 1 )
#     else
        CALL lbc( zbetup(1,1,1), 1, 1, 1, 1, jpk, 1 )
        CALL lbc( zbetdo(1,1,1), 1, 1, 1, 1, jpk, 1 )
#     endif

C     
C
C 3. monotonic flux in the i direction, i.e. paa
C ----------------------------------------------
C
      DO jk=1,jpkm1
        DO jj=2,jpjm1
          DO ji=2,jpim1
            zc= paa(ji,jj,jk)
            IF (zc.ge.0.) THEN
                za= MIN(1.D0,zbetdo(ji,jj,jk),zbetup(ji+1,jj,jk))
                paa(ji,jj,jk) = za * zc
              ELSE
                zb= MIN(1.D0,zbetup(ji,jj,jk),zbetdo(ji+1,jj,jk))
                paa(ji,jj,jk) = zb * zc
            ENDIF
          END DO
        END DO
      END DO
C
Cek   CALL bord3(paa)
#     ifdef key_mpp
        CALL mpplnk( paa(1,1,1), 1, 1 )
#     else
        CALL lbc( paa(1,1,1), 1, 1, 1, 1, jpk, 1 )
#     endif

C
C
C 4. monotonic flux in the j direction, i.e. pbb
C ----------------------------------------------
C
      DO jk=1,jpkm1
        DO jj=2,jpjm1
          DO ji=2,jpim1
            zc= pbb(ji,jj,jk)
            IF (zc.ge.0.) THEN
                za= MIN(1.D0,zbetdo(ji,jj,jk),zbetup(ji,jj+1,jk))
                pbb(ji,jj,jk) = za * zc
              ELSE
                zb= MIN(1.D0,zbetup(ji,jj,jk),zbetdo(ji,jj+1,jk))
                pbb(ji,jj,jk) = zb * zc
            ENDIF
          END DO
        END DO
      END DO
C
Cek   CALL bord3(pbb)
#     ifdef key_mpp
        CALL mpplnk( pbb(1,1,1), 1, 1 )
#     else
        CALL lbc( pbb(1,1,1), 1, 1, 1, 1, jpk, 1 )
#     endif

C
C
C 5. monotonic flux in the k direction, i.e. pcc
C ----------------------------------------------
C
      DO jk=2,jpkm1
        DO jj=2,jpjm1
          DO ji=2,jpim1
            zc= pcc(ji,jj,jk)
            IF (zc.ge.0.) THEN
                za= MIN(1.D0,zbetdo(ji,jj,jk),zbetup(ji,jj,jk-1))
                pcc(ji,jj,jk) = za * zc
              ELSE
                zb= MIN(1.D0,zbetup(ji,jj,jk),zbetdo(ji,jj,jk-1))
                pcc(ji,jj,jk) = zb * zc
            ENDIF
          END DO
        END DO
      END DO
C
Cek   CALL bord3(pcc)
#     ifdef key_mpp
        CALL mpplnk( pcc(1,1,1), 1, 1 )
#     else
        CALL lbc( pcc(1,1,1), 1, 1, 1, 1, jpk, 1 )
#     endif

C
#endif
C
C
      RETURN
      END
