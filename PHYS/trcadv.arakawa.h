CCC $Id: trcadv.arakawa.h,v 1.2 2009-09-11 09:20:56 cvsogs01 Exp $
CCC
CCC                      trcadv.arakawa.h
CCC                     ******************
CCC
CC   defined key : 'key_trc_arakawa'
CC   ============
CC
CC  Purpose :
CC  ---------
CC     Compute the now trend due to the horizontal and vertical
CC      advection of tracers 
CC     (tr) and add it to the general trend of passive tracer equations.
CC
CC
CC   Method :
CC   -------
CC      This routine compute not exactly the advection but the
CC      transport term, i.e.  div(utr).
CC
CC      The now horizontal advection of tracers is given by:
CC              dx(trn un) + dy(trn vn) = 1/(e1t*e2t*e3t)
CC                                           ( di-1( e2u*e3u mi(trn) )
CC                                           + dj-1( e1v*e3v mj(trn) )
CC         N.B. if neither 'key_s_coord' activated (z-coord.) a
CC      simplification by e3 is done in the previous expression.
CC
CC      Add this trend now to the general trend of tracer (ta,sa):
CC		tra = tra + dx(trn un) + dy(trn vn)
CC
CC	The now vertical advection of tracers is given by:
CC            dz(trn wn) = 1/( e1t*e2t*e3t) dk+1( e1t e2t wn mk(trn) )
CC         since horizontal scale factor do not depend on the level
CC         these trends reduced to:
CC            dz(trn wn) = 1/e3t dk+1( wn mk(trn) )
CC
CC	Add this trend now to the general trend of tracer (tra):
CC            tra = tra + dz(trn wn)
CC
CC      'key_trc_diatrd' activated: the now horizontal advection trend
CC	is saved for futher diagnostics.
CC
CC      macro-tasked on tracer slab (jn-loop)
CC
CC
CC   Input :
CC   -----
CC      argument
CC              ktask           : task identificator
CC              kt              : time step
CC      common
CC            /comcoh/          : orthogonal curvilinear coordinates
CC                                and scale factors
CC            /cottrp/ trn,tra  : present fields (now) and next fields
CC                                (after) for passive tracer
CC            /comtsk/          : multitasking
CC
CC   Output :
CC   ------
CC      common
CC            /cottrp/ tra      : general tracer trend increased by the
CC                                now horizontal tracer advection trend
CC            /cottrdD/ trtrd     : now horizontal tracer advection trend
CC                                (if 'key_trc_diatrd' is activated)
CC
CC   Modifications:
CC   --------------
CC      original : 87-06 (pa-dl)
CC      addition : 91-11 (G. Madec)
CC      addition : 95-11 (G. Madec) suppress volumetric scale factors
CC      addition : 96-01 (G. Madec) terrain following coordinates
CC                                  statement function for e3
CC                                  suppression of common work arrays
CC               : 97-04 (M.A. Foujols) adapted for passive tracer
CC                                  suppresion of ice specification 
CC      modifications : 00-05 (MA Foujols) add lbc for tracer trends
CC      modifications : 00-10 (MA Foujols and E. Kestenare) INCLUDE
CC                             instead of routine and had and zad concatenation
CC----------------------------------------------------------------------
      INTEGER  ji,jj,jk,jn
      REAL(8) zbtr,ztra,zbe1ru,zbe2rv,zhw,ze3tr
      REAL(8) zwx(jpi,jpj),zwy(jpi,jpj),zwz(jpi,jpk)
CC----------------------------------------------------------------------
CC statement functions
CC ===================

!               #include "stafun.h"

CCC---------------------------------------------------------------------
CCC  OPA8.1, LODYC (1998)
CCC---------------------------------------------------------------------
C
C
C Tracer slab
C ===========
C
      DO jn=1,jptra
C
C Horizontal slab
C ===============
C
        DO jk=1,jpkm1
C
C
C 1. Horizontal advection of passive tracers
C ------------------------------------------
C
C Tracer flux at u and v-points
C
          DO jj=1,jpjm1
            DO ji=1,jpim1
#if defined key_s_coord
              zbe1ru = 0.5 * e2u(ji,jj)*fse3u(ji,jj,jk) * un(ji,jj,jk)
              zbe2rv = 0.5 * e1v(ji,jj)*fse3v(ji,jj,jk) * vn(ji,jj,jk)
#  else
              zbe1ru     = 0.5 * e2u(ji,jj) * un(ji,jj,jk)
              zbe2rv     = 0.5 * e1v(ji,jj) * vn(ji,jj,jk)
#endif
              zwx(ji,jj) = zbe1ru *
     $            ( trn(ji,jj,jk,jn) + trn(ji+1,jj  ,jk,jn) )
              zwy(ji,jj) = zbe2rv *
     $            ( trn(ji,jj,jk,jn) + trn(ji  ,jj+1,jk,jn) )
            END DO  
          END DO  
C
C Tracer flux divergence at t-point added to the general trend
C
          DO jj=2,jpjm1
            DO ji=2,jpim1
#if defined key_s_coord
              zbtr= 1. / ( e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk) )
#  else
              zbtr= 1. / ( e1t(ji,jj)*e2t(ji,jj) )
#endif
              ztra = - zbtr * (  ( zwx(ji,jj) - zwx(ji-1,jj  ) )
     $            + ( zwy(ji,jj) - zwy(ji  ,jj-1) )  )
              tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn) + ztra
#if defined key_trc_diatrd
              trtrd(ji,jj,jk,jn,1) = - zbtr * 
     $                             ( zwx(ji,jj) - zwx(ji-1,jj) )
              trtrd(ji,jj,jk,jn,2) = - zbtr * 
     $                             ( zwy(ji,jj) - zwy(ji,jj-1) )
#endif
            END DO  
          END DO  
C
        END DO 

C Lateral boundary conditions on trtrd:
#      if defined key_trc_diatrd
#         ifdef key_mpp
        CALL mpplnk( trtrd(1,1,1,jn,1), 1, 1 )
        CALL mpplnk( trtrd(1,1,1,jn,2), 1, 1 )
#         else      
        CALL lbc( trtrd(1,1,1,jn,1), 1, 1, 1, 1, jpk, 1 )
        CALL lbc( trtrd(1,1,1,jn,2), 1, 1, 1, 1, jpk, 1 )
#         endif
#      endif
C
C Vertical slab
C =============
C
        DO jj=1,jpj
C
C
C 1. Vertical advection of tracers
C --------------------------------
          zwz=0.
C
C 1.1 Tracer flux at w-point
C
          DO jk=2,jpk
            DO ji=2,jpim1
              zhw = 0.5 * wn(ji,jj,jk)
              zwz(ji,jk) = zhw*( trn(ji,jj,jk,jn) + trn(ji,jj,jk-1,jn) )
            END DO  
          END DO  
C
C 1.2 Surface value: zero 
C
          DO ji=2,jpim1
            zwz(ji,1) = 0.e0
          END DO  
C
C
C 1.3 Tracer flux divergence at t-point added to the general trend
C
          DO jk=1,jpkm1
            DO ji=2,jpim1
              ze3tr = 1./fse3t(ji,jj,jk)
              ztra = - ze3tr * ( zwz(ji,jk) - zwz(ji,jk+1) )
              tra(ji,jj,jk,jn) =  tra(ji,jj,jk,jn) + ztra
#if defined key_trc_diatrd
              trtrd(ji,jj,jk,jn,3) = ztra
#endif
            END DO  
          END DO  
C
        END DO 
#      if defined key_trc_diatrd
#         ifdef key_mpp
        CALL mpplnk( trtrd(1,1,1,jn,3), 1, 1 )
#         else      
        CALL lbc( trtrd(1,1,1,jn,3), 1, 1, 1, 1, jpk, 1 )
#         endif
#      endif
C
C End of tracer slab
C ==================
C
       END DO

