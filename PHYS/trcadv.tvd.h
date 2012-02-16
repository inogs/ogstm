CCC $Id: trcadv.tvd.h,v 1.2 2009-09-11 09:20:56 cvsogs01 Exp $
CCC
CCC                      trcadv.tvd.h
CCC                     **************
CCC
CC   defined key : 'key_trc_tvd'
CC   ============
CC
CC  PURPOSE :
CC ---------
CC     Compute the now trend due to total advection of tracers
CC     and add it to the general trend of tracer equations:
CC
CC   METHOD :
CC   -------
CC      this ROUTINE compute not exactly the advection but the
CC      transport term, i.e.  div(u*tra).
CC
CC   centered sheme withe corrected flux (monotonic corection)
CC   *********************************************************
CC
CC      cf REFERENCE
CC
CC      note: - this advection scheme needs a leap-frog time scheme
CC
CC   remarks :
CC   -------  
CC      IF 'key_trc_diatrd' key is activated, the
CC      advection trend is saved for futher diagnostics.
CC
CC      multitasked on tracer (jn-loop)
CC
CC      COMMON
CC            /comcoo/          : orthogonal curvilinear coordinates
CC                                and scale factors
CC            /cottrp/          : present and next fields for passive
CC                              : tracers
CC            /comtsk/          : multitasking
CC
CC   OUTPUT :
CC   ------
CC      COMMON
CC            /cottrp/ tra      : general tracer trend increased by the
CC                                now horizontal tracer advection trend
CC            /cottrd/ trtrd    : now horizontal tracer advection trend
CC                                (IF 'key_trc_diatrd' key is activated)
CC
CC
CC   EXTERNAL :                  nonosc, mpplnk, lbc
CC   --------
CC

      USE TIME_MANAGER
      INTEGER ji,jj,jk,jn
      REAL(8) zti(jpi,jpj,jpk)
      REAL(8) zkx(jpi,jpj,jpk),zky(jpi,jpj,jpk),zkz(jpi,jpj,jpk)
      REAL(8) z2dtt,zbt,zeu,zev,zew,z2
C
CC----------------------------------------------------------------------
CC statement functions
CC ===================

#include "stafun.h"

      IF ((neuler.eq.0).and.(kt.eq.TimeStepStart)) THEN
           z2=1.
        ELSE
           z2=2.
      ENDIF

C
C tracer loop parallelized (macrotasking)
C =======================================
C
      DO 1000 jn = ktask,jptra,ncpu
C
C 1. initialization
C -----------------
C
        zkx = 0.
        zky = 0.
        zkz = 0.
        zti = 0.
C
C
C 2. upstream advection with initial mass fluxes and intermediate update 
C ----------------------------------------------------------------------
C
C 2.1 upstream tracer flux in the i and j direction
C
        DO jk=1,jpk
          DO jj = 1,jpjm1
            DO ji = 1,jpim1
              zeu = e2u(ji,jj)* fse3u(ji,jj,jk)* un(ji,jj,jk)
              zev = e1v(ji,jj)* fse3v(ji,jj,jk)* vn(ji,jj,jk)
              if (un(ji,jj,jk).gt.0.) then 
                  zkx(ji,jj,jk) = trb(ji,jj,jk,jn)   * zeu
               else 
                  zkx(ji,jj,jk) = trb(ji+1,jj,jk,jn) * zeu
              endif 
              if (vn(ji,jj,jk).gt.0.) then
                  zky(ji,jj,jk) = trb(ji,jj,jk,jn)   * zev
               else 
                  zky(ji,jj,jk) = trb(ji,jj+1,jk,jn) * zev
              endif 
            END DO
          END DO
        END DO
C
C 2.2 upstream tracer flux in the k direction
C
        DO jk = 2,jpk
          DO jj = 1,jpj
            DO ji = 1,jpi
              zew = e1t(ji,jj)*   e2t(ji,jj)*    wn(ji,jj,jk)
              if(wn(ji,jj,jk).gt.0.) then
                 zkz(ji,jj,jk) = trb(ji,jj,jk,jn) * zew
               else
                 zkz(ji,jj,jk) = trb(ji,jj,jk-1,jn) * zew
              endif 
            END DO
          END DO
        END DO

C
C 2.3 total advective trend
C
        DO jk=1,jpkm1
          DO jj=2,jpjm1
            DO ji=2,jpim1
              zbt = e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk)
              zti(ji,jj,jk)= - ( zkx(ji,jj,jk) - zkx(ji-1,jj,jk)
     $                          +zky(ji,jj,jk) - zky(ji,jj-1,jk)
     $                          +zkz(ji,jj,jk) - zkz(ji,jj,jk+1)
     $                          ) / zbt
#if defined key_trc_diatrd
              trtrd(ji,jj,jk,jn,1) =
     $              - ( zkx(ji,jj,jk) - zkx(ji-1,jj,jk) ) / zbt
              trtrd(ji,jj,jk,jn,2) =
     $              - ( zky(ji,jj,jk) - zky(ji,jj-1,jk) ) / zbt
              trtrd(ji,jj,jk,jn,3) =
     $              - ( zkz(ji,jj,jk) - zkz(ji,jj,jk+1) ) / zbt
#endif
            END DO
          END DO
        END DO

C
C 2.4 update and guess with monotonic sheme
C
        DO jk=1,jpkm1
          z2dtt = z2 * rdt
          DO jj=2,jpjm1
            DO ji=2,jpim1
              tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn) + zti(ji,jj,jk)
              zti(ji,jj,jk) = (trb(ji,jj,jk,jn)+z2dtt*zti(ji,jj,jk))
     $                          *tmask(ji,jj,jk)
            END DO
          END DO
        END DO
C
C
C
C ... Lateral boundary conditions on zti
#ifdef key_mpp
C
C   ... Mpp : export boundary values to neighboring processors
C
        CALL mpplnk( zti, 1, 1 )
#else
C
C   ... T-point, 3D array, full local arrays zti are initialised
C       
C
        CALL lbc( zti, 1, 1, 1, 1, jpk, 1 )
C
#endif
C
C
C 3. antidiffusive flux : high order minus low order
C --------------------------------------------------
C
C 3.1 antidiffusive flux on i and j
C
        DO jk=1,jpk
          DO jj=1,jpjm1
            DO ji=1,jpim1
              zkx(ji,jj,jk) = e2u(ji,jj)*fse3u(ji,jj,jk)*un(ji,jj,jk)
     $              *(trn(ji,jj,jk,jn)+trn(ji+1,jj,jk,jn)) / 2.
     $                  - zkx(ji,jj,jk)
              zky(ji,jj,jk) = e1v(ji,jj)*fse3v(ji,jj,jk)*vn(ji,jj,jk)
     $              *(trn(ji,jj,jk,jn)+trn(ji,jj+1,jk,jn)) / 2.
     $                  - zky(ji,jj,jk)
            END DO
          END DO
        END DO
C
C 3.2 antidiffusive flux on k
C
        DO jk=2,jpk
          DO jj=1,jpj
            DO ji=1,jpi
              zkz(ji,jj,jk) = e1t(ji,jj)*e2t(ji,jj) * wn(ji,jj,jk)
     $              *(trn(ji,jj,jk,jn)+trn(ji,jj,jk-1,jn)) / 2.
     $                  - zkz(ji,jj,jk)
            END DO
          END DO
        END DO
C
C 3.3 boundary conditions
C
        DO jj=1,jpj
          DO ji=1,jpi
            zkz(ji,jj,1) = 0.
          END DO
        END DO
C
#ifdef key_mpp
C
C   ... Mpp : export boundary values to neighboring processors
C
          CALL mpplnk( zkx, 1, 1 )
          CALL mpplnk( zky, 1, 1 )
          CALL mpplnk( zkz, 1, 1 )
#else
C
C   ... T-point, 3D array, full local array zk[xy] are initialised
C       
C
          CALL lbc( zkx, 1, 1, 1, 1, jpk, 1 )
          CALL lbc( zky, 1, 1, 1, 1, jpk, 1 )
          CALL lbc( zkz, 1, 1, 1, 1, jpk, 1 )
C
#endif
C
C
C 4. monotonicity algorithm
C -------------------------
C
        z2dtt = z2*rdt
        CALL nonosc(trb(1,1,1,jn),zkx,zky,zkz,zti,z2dtt)
C
C
C 5. final trend with corrected fluxes
C ------------------------------------
C
        DO jk=1,jpkm1
          DO jj=2,jpjm1
            DO ji=2,jpim1
              zbt = e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk)
              tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn) - (
     $                          zkx(ji,jj,jk) - zkx(ji-1,jj,jk)
     $                         +zky(ji,jj,jk) - zky(ji,jj-1,jk)
     $                         +zkz(ji,jj,jk) - zkz(ji,jj,jk+1)
     $                          )  / zbt
#if defined key_trc_diatrd
              trtrd(ji,jj,jk,jn,1) = trtrd(ji,jj,jk,jn,1) -
     $               ( zkx(ji,jj,jk) - zkx(ji-1,jj,jk) ) / zbt
              trtrd(ji,jj,jk,jn,2) = trtrd(ji,jj,jk,jn,2) -
     $               ( zky(ji,jj,jk) - zky(ji,jj-1,jk) ) / zbt
              trtrd(ji,jj,jk,jn,3) = trtrd(ji,jj,jk,jn,3) -
     $               ( zkz(ji,jj,jk) - zkz(ji,jj,jk+1) ) / zbt
#endif
            END DO
          END DO
        END DO

#      if defined key_trc_diatrd
#         ifdef key_mpp
        CALL mpplnk( trtrd(1,1,1,jn,1), 1, 1 )
        CALL mpplnk( trtrd(1,1,1,jn,2), 1, 1 )
        CALL mpplnk( trtrd(1,1,1,jn,3), 1, 1 )
#         else      
        CALL lbc( trtrd(1,1,1,jn,1), 1, 1, 1, 1, jpk, 1 )
        CALL lbc( trtrd(1,1,1,jn,2), 1, 1, 1, 1, jpk, 1 )
        CALL lbc( trtrd(1,1,1,jn,3), 1, 1, 1, 1, jpk, 1 )
#         endif
#      endif

C
C END of tracer loop
C ==================
C
 1000 CONTINUE
C
