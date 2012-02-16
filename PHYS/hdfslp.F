CC $Header: /cvsroot/opatm-bfm/opa_model/OPA/hdfslp.F,v 1.2 2009-09-11 09:20:56 cvsogs01 Exp $

      SUBROUTINE hdfslp ( ktask, kt )
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE hdfslp
CCC                     ******************
CCC
CCC  Purpose :
CCC  --------
CCC     Compute the slopes of neutral surface (slope of isopycnal sur-
CCC	faces referenced locally) ('key_trahdfiso' defined).
CCC
CC   Method :
CC   -------
CC      The slope in the i-direction is computed at u- and w-points
CC	(uslp, wslpi) and the slope in the j-direction is computed at 
CC	v- and w-points (vslp, wslpj).
CC	They are bounded by 1/100 over the whole ocean, and within the
CC	surface layer they are bounded by the distance to the surface
CC	( slope<= depth/l  where l is the length scale of horizontal
CC	diffusion (here, aht=2000m2/s ==> l=20km with a typical velocity
CC	of 10cm/s)
CC
CC      A horizontal shapiro filter is applied to the slopes
CC
CC	'key_s_coord' defined: add to the previously computed slopes
CC	the slope of the model level surface.
CC
CC      macro-tasked on horizontal slab (jk-loop)  (2, jpk-1)
CC	[slopes already set to zero at level 1, and to zero or the ocean
CC	bottom slope ('key_s_coord' defined) at level jpk in inihdf]
CC
CC   Input :
CC   ------
CC      argument
CC              ktask           : task identificator
CC              kt              : time step
CC      common
CC            /comcoh/          : scale factors
CC            /comask/          : masks, bathymetry
CC            /comnow/          : now fields 
CC            /comtsk/          : multitasking
CC
CC   Output :
CC   -------
CC      common
CC            /comiso/uslp(),   : i- and j-slopes of neutral surfaces
CC                    vslp()      at u- and v-points, resp.
CC		      wslpi(),  : i- and j-slopes of neutral surfaces
CC		      wslpj()     at w-points. 
CC
CC   External :       mpplnk, mpplnk2  or  lbc

       USE myalloc
        IMPLICIT NONE


CC----------------------------------------------------------------------
CC local declarations

      INTEGER ktask, kt
#if defined key_trahdfiso
      INTEGER ji, jj, jk
C
      REAL(8) zwy(jpi,jpj), zwz(jpi,jpj), zww(jpi,jpj)
      REAL(8) zeps, zmg, zm05g, zcoef1, zcoef2
      REAL(8) zau, zbu, zav, zbv, zslopu, zslopv
      REAL(8) zai, zbi, zaj, zbj, zslopm
      REAL(8) zcofu, zcofv, zcofw, z1u, z1v, z1w, zclmaxw, zclmaxz
      REAL(8) zcof, zfr, zfri, zfrj, zfrk
CC----------------------------------------------------------------------
CC statement functions
CC ===================

#include "stafun.h"

CCC---------------------------------------------------------------------
CCC  OPA8, LODYC (1997)
CCC---------------------------------------------------------------------
C
C
C 0. Local constant initialization
C --------------------------------
C
      zeps  =  1.e-20
      zmg   = -1.0 / g
      zm05g = -0.5 / g
C
C
C Vertical slab
C =============
C
      DO jk = ktask+1, jpkm1
C
        DO jj = 1, jpj
          DO ji = 1, jpi
            zwz(ji,jj) = 0.0
            zww(ji,jj) = 0.0
          END DO
        END DO
C
C I.  slopes at u and v point
C ===========================
C
C
C I.1. Slopes of isopycnal surfaces
C ---------------------------------
C uslp = d/di( rdn ) / d/dz( rdn )
C vslp = d/dj( rdn ) / d/dz( rdn )
C
C ... Local vertical density gradient evaluated from N^2
C     zwy = d/dz(rdn)= - ( rdn ) / g * mk(bn2n) -- at t point
C
        DO jj = 1, jpj
          DO ji = 1, jpi
            zwy(ji,jj) = zmg * ( rdn(ji,jj,jk) + 1. )
     $                 * ( bn2n(ji,jj,jk) + bn2n(ji,jj,jk+1) )
     $                 / max( tmask(ji,jj,jk) + tmask (ji,jj,jk+1), 1. )
          END DO
        END DO
C
C ... Slope at u and v points
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
C   ... horizontal and vertical density gradient at u- and v-points
            zau = 1. / e1u(ji,jj) * ( rdn(ji+1,jj,jk) - rdn(ji,jj,jk) )
            zav = 1. / e2v(ji,jj) * ( rdn(ji,jj+1,jk) - rdn(ji,jj,jk) )
            zbu = 0.5*( zwy(ji,jj)+zwy(ji+1,jj) )
            zbv = 0.5*( zwy(ji,jj)+zwy(ji,jj+1) )
C   ... bound the slopes: abs(zw.)<= 1/100 and zb..<0
C                         kxz max= ah slope max =< e1 e3 /(pi**2 2 dt)
            zbu = min( zbu, -100.*abs(zau),
     $                   -7.e+3/fse3u(ji,jj,jk)*abs(zau) )
            zbv = min( zbv, -100.*abs(zav),
     $                   -7.e+3/fse3v(ji,jj,jk)*abs(zav) )
C   ... uslp and vslp output in zwz and zww, resp.
            zwz (ji,jj) = zau / ( zbu - zeps ) * umask (ji,jj,jk)
            zww (ji,jj) = zav / ( zbv - zeps ) * vmask (ji,jj,jk)
          END DO
        END DO
C
C ... lateral boundary conditions on zww and zwz
#ifdef key_mpp
C   ... Mpp: export boundary values to neighbouring processors
CCC 25/06/2003 Cristiano - It has been adjusted (add a parameter), but I am not sure of the value
CCC        CALL mpplnk2( zwz, 2 )
CCC        CALL mpplnk2( zww, 3 )
        CALL mpplnk2( zwz, 2, 1)
        CALL mpplnk2( zww, 3, 1)

#  else
C   ... mono- or macro-tasking: U- and V-points, 2D array, no slab
            CALL lbc2( zwz, 2, 0, 1, 1, 1, 1 )
            CALL lbc2( zww, 3, 0, 1, 1, 1, 1 )
#endif
C
C ... Shapiro filter applied in the horizontal direction
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zcofu = 1. / 16.
            zcofv = 1. / 16.
C   ... uslop
            uslp(ji,jj,jk) = zcofu * ( zwz(ji-1,jj-1) + zwz(ji+1,jj-1)
     $                               + zwz(ji-1,jj+1) + zwz(ji+1,jj+1)
     $                          + 2.*( zwz(ji  ,jj-1) + zwz(ji-1,jj  )
     $                                +zwz(ji+1,jj  ) + zwz(ji  ,jj+1) )
     $                            + 4.*zwz(ji  ,jj  )   )
C   ... vslop
            vslp(ji,jj,jk) = zcofv * ( zww(ji-1,jj-1) + zww(ji+1,jj-1)
     $                               + zww(ji-1,jj+1) + zww(ji+1,jj+1)
     $                          + 2.*( zww(ji  ,jj-1) + zww(ji-1,jj  )
     $                                +zww(ji+1,jj  ) + zww(ji  ,jj+1) )
     $                            + 4.*zww(ji,jj  )      )
          END DO
        END DO
C
C ... decrease along coastal boundaries
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            z1u = ( 1. - abs( fmask(ji,jj,jk) - fmask(ji-1,jj,jk) )*.5 )
            z1v = ( 1. - abs( fmask(ji,jj,jk) - fmask(ji,jj-1,jk) )*.5 )
            z1w =           ( tmask(ji,jj,jk) + tmask(ji,jj,jk+1) )*.5
            uslp(ji,jj,jk) = uslp(ji,jj,jk)       * z1v * z1w
            vslp(ji,jj,jk) = vslp(ji,jj,jk) * z1u       * z1w
          END DO
        END DO
C
#ifdef key_ice
C
C ... zero slope at the vicinity of freezing point
        IF (jk.LE.jpk) THEN
            DO jj = 1, jpjm1
              DO ji = 1, jpim1
                zcof = fzptb(ji,jj) + 0.1
                zfr  = tmask(ji,jj,jk)
     $               * max( 0., sign( 1., zcof - tb(ji,jj,jk) ) )
                zfri = tmask(ji+1,jj,jk)
     $               * max( 0., sign( 1., zcof - tb(ji+1,jj,jk) ) )
                zfrj = tmask(ji,jj+1,jk)
     $               * max( 0., sign( 1., zcof - tb(ji,jj+1,jk) ) )
                z1u = ( 1. - zfr ) * ( 1. - zfri )
                z1v = ( 1. - zfr ) * ( 1. - zfrj )
                uslp(ji,jj,jk) = uslp(ji,jj,jk) * z1u
                vslp(ji,jj,jk) = vslp(ji,jj,jk) * z1v
              END DO
            END DO
        ENDIF
#endif
C
#if defined key_s_coord
C
C I.2 Add the slope of level surfaces
C -----------------------------------
C 'key_s_coord' defined but not 'key_trahdfiso' the computation is done
C in inihdf, hdfslp never called
C 'key_s_coord' and 'key_trahdfiso' defined, the slope of level surfaces
C is added to the slope of isopycnal surfaces.
C c a u t i o n : minus sign as fsdep has positive value 
C
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            uslp(ji,jj,jk) = uslp(ji,jj,jk) - 1. / e1u(ji,jj)
     $                     * ( fsdept(ji+1,jj,jk) - fsdept(ji,jj,jk) )
            vslp(ji,jj,jk) = vslp(ji,jj,jk) - 1. / e2v(ji,jj)
     $                     * ( fsdept(ji,jj+1,jk) - fsdept(ji,jj,jk) )
          END DO
        END DO
C
#endif
C
C I.3 Surface boundary layer (bound the physical slope)
C --------------------------
C maximum slope within the surface layer is bounded by the distance
C to the surface ( slopemax<= depth/l  where l is the length scale
C off horizontal diffusion (here, aht=2000m2/s ==> l=20km with
C typical velocity of 10cm/s)
C
        zslopu = fsdept(1,2,jk) / 2.e+4
        zslopv = zslopu
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
#    if defined key_s_coord
            zslopu = 0.5*( fsdept(ji,jj,jk)+fsdept(ji+1,jj,jk) ) / 2.e+4
            zslopv = 0.5*( fsdept(ji,jj,jk)+fsdept(ji,jj+1,jk) ) / 2.e+4
#    endif
            zclmaxw        = max( uslp(ji,jj,jk), -zslopu )
            uslp(ji,jj,jk) = min( zclmaxw       ,  zslopu )
            zclmaxz        = max( vslp(ji,jj,jk), -zslopv )
            vslp(ji,jj,jk) = min( zclmaxz       ,  zslopv )
          END DO
        END DO
C
C
C
C II. Computation of slopes at w point
C ====================================
C
C
C II.1 Slopes of isopycnal surfaces
C ---------------------------------
C wslpi = mij( d/di( rdn ) / d/dz( rdn )
C wslpj = mij( d/dj( rdn ) / d/dz( rdn )
C
C
C ... Local vertical density gradient evaluated from N^2
C     zwy = d/dz(rdn)= - mk ( rdn ) / g * bn2n -- at w point
        DO jj = 1, jpj
          DO ji = 1, jpi
            zwy (ji,jj) = zm05g * bn2n (ji,jj,jk) *
     $          ( rdn (ji,jj,jk) + rdn (ji,jj,jk-1) + 2. )
          END DO
        END DO
C
C ... Slope at w point
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
C   ... horizontal density i-gradient at w-points
            zcoef1 = max( zeps, umask(ji-1,jj,jk  )+umask(ji,jj,jk  )
     $                         +umask(ji-1,jj,jk-1)+umask(ji,jj,jk-1) )
            zcoef1 = 1. / ( zcoef1 * e1t (ji,jj) )
            zai = zcoef1 * (
     $   ( rdn(ji+1,jj, jk ) - rdn( ji ,jj, jk ) ) * umask( ji ,jj, jk )
     $ + ( rdn(ji+1,jj,jk-1) - rdn( ji ,jj,jk-1) ) * umask( ji ,jj,jk-1)
     $ + ( rdn( ji ,jj,jk-1) - rdn(ji-1,jj,jk-1) ) * umask(ji-1,jj,jk-1)
     $ + ( rdn( ji ,jj, jk ) - rdn(ji-1,jj, jk ) ) * umask(ji-1,jj, jk )
     $                      ) * tmask (ji,jj,jk)
C   ... horizontal density j-gradient at w-points
            zcoef2 = max( zeps, vmask(ji,jj-1,jk  )+vmask(ji,jj,jk-1)
     $                         +vmask(ji,jj-1,jk-1)+vmask(ji,jj,jk  ) )
            zcoef2 = 1.0 / ( zcoef2 *  e2t (ji,jj) )
            zaj = zcoef2 * (
     $   ( rdn(ji,jj+1, jk ) - rdn(ji, jj , jk ) ) * vmask(ji, jj , jk )
     $ + ( rdn(ji,jj+1,jk-1) - rdn(ji, jj ,jk-1) ) * vmask(ji, jj ,jk-1)
     $ + ( rdn(ji, jj ,jk-1) - rdn(ji,jj-1,jk-1) ) * vmask(ji,jj-1,jk-1)
     $ + ( rdn(ji ,jj , jk ) - rdn(ji,jj-1, jk ) ) * vmask(ji,jj-1, jk )
     $                      ) * tmask (ji,jj,jk)
C   ... bound the slopes: abs(zw.)<= 1/100 and zb..<0.
C                         static instability:
C                         kxz max= ah slope max =< e1 e3 /(pi**2 2 dt)
            zbi = min ( zwy (ji,jj),- 100.*abs(zai),
     $                   -7.e+3/fse3w(ji,jj,jk)*abs(zai) )
            zbj = min ( zwy (ji,jj), -100.*abs(zaj),
     $                   -7.e+3/fse3w(ji,jj,jk)*abs(zaj) )
C   ... wslpi and wslpj output in zwz and zww, resp.
            zwz (ji,jj) = zai / ( zbi - zeps) * tmask (ji,jj,jk)
            zww (ji,jj) = zaj / ( zbj - zeps) * tmask (ji,jj,jk)
          END DO
        END DO
C
C ... lateral boundary conditions on zwz and zww
#ifdef key_mpp
C   ... Mpp: export boundary values to neighbouring processors
CCC 25/06/2003 Cristiano - It has been adjusted (add a parameter), but I am not sure of the value
CCC        CALL mpplnk2( zwz, 1 )
CCC        CALL mpplnk2( zww, 1 )
        CALL mpplnk2( zwz, 1, 1)
        CALL mpplnk2( zww, 1, 1)
#  else
C   ... mono- or macro-tasking: T-point, 2D array, no slab
            CALL lbc2( zwz, 1, 0, 1, 1, 1, 1 )
            CALL lbc2( zww, 1, 0, 1, 1, 1, 1 )
#endif
C
C ... Shapiro filter applied in the horizontal direction
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zcofw = tmask(ji,jj,jk)/16.
            wslpi(ji,jj,jk) = (       zwz(ji-1,jj-1) + zwz(ji+1,jj-1)
     $                         +      zwz(ji-1,jj+1) + zwz(ji+1,jj+1)
     $                         + 2.*( zwz(ji  ,jj-1) + zwz(ji-1,jj  )
     $                               +zwz(ji+1,jj  ) + zwz(ji  ,jj+1) )
     $                         + 4.*  zwz(ji,jj  )       ) * zcofw
            wslpj(ji,jj,jk) = (       zww(ji-1,jj-1) + zww(ji+1,jj-1)
     $                         +      zww(ji-1,jj+1) + zww(ji+1,jj+1)
     $                         + 2.*( zww(ji  ,jj-1) + zww(ji-1,jj  )
     $                               +zww(ji+1,jj  ) + zww(ji  ,jj+1) )
     $                         + 4.*  zww(ji,jj  )       ) * zcofw
          END DO
        END DO
C
C ... decrease the slope along the boundaries
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            z1u = ( umask(ji,jj,jk) + umask(ji-1,jj,jk) ) *.5
            z1v = ( vmask(ji,jj,jk) + vmask(ji,jj-1,jk) ) *.5
            wslpi(ji,jj,jk) = wslpi(ji,jj,jk) * z1u * z1v
            wslpj(ji,jj,jk) = wslpj(ji,jj,jk) * z1u * z1v
          END DO
        END DO
C
C ... bottom value set to 0.
        IF (jk.EQ.jpk) THEN
            DO jj = 1, jpj
              DO ji = 1, jpi
                wslpi(ji,jj,jpk) = 0.e0
                wslpj(ji,jj,jpk) = 0.e0
              END DO
            END DO
        ENDIF
C
#ifdef key_ice
C
C ... zero slope in the vicinity of freezing point
        IF (jk.LE.jpkm1) THEN
            DO jj = 1, jpj
              DO ji = 1, jpi
                zcof = fzptb(ji,jj) + 0.1
                zfr  = tmask(ji,jj,jk)
     $               * max( 0., sign( 1., zcof - tb(ji,jj,jk  ) ) )
                zfrk = tmask(ji,jj,jk-1)
     $               * max( 0., sign( 1., zcof - tb(ji,jj,jk-1) ) )
                z1w = ( 1. - zfr ) * ( 1. - zfrk )
                wslpi(ji,jj,jk) = wslpi(ji,jj,jk) * z1w
                wslpj(ji,jj,jk) = wslpj(ji,jj,jk) * z1w
              END DO
            END DO
        ENDIF
#endif
C
#if defined key_s_coord
C
C II.2 Slope of level surfaces
C ----------------------------
C 'key_s_coord' defined but not 'key_trahdfiso' the computation is done
C in inihdf, hdfslp never called
C 'key_s_coord' and 'key_trahdfiso' defined, the slope of level surfaces
C is added to the slope of isopycnal surfaces.
C
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            wslpi(ji,jj,jk) = wslpi(ji,jj,jk) - 1. / e1t(ji,jj)
     $                     * ( fsdepuw(ji+1,jj,jk) - fsdepuw(ji,jj,jk) )
            wslpj(ji,jj,jk) = wslpj(ji,jj,jk) - 1. / e2t(ji,jj)
     $                     * ( fsdepvw(ji,jj+1,jk) - fsdepvw(ji,jj,jk) )
          END DO
        END DO
C
#endif
C
C II.3. Surface boundary layer (bound the physical slope)
C ----------------------------
C maximum slope within the surface layer is bounded by the distance
C to the surface ( slopemax<= depth/l  where l is the length scale
C off horizontal diffusion (here, aht=2000m2/s ==> l=20km with
C typical velocity of 10cm/s)
C
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zslopm = fsdepw(ji,jj,jk) / 2.e+4
            zclmaxw         = max( wslpi(ji,jj,jk), -zslopm )
            wslpi(ji,jj,jk) = min( zclmaxw        ,  zslopm )
            zclmaxz         = max( wslpj(ji,jj,jk), -zslopm )
            wslpj(ji,jj,jk) = min( zclmaxz        ,  zslopm )
          END DO
        END DO
C
C
C End of task
C ===========
C
       END DO
C
C
C III Lateral boundary conditions on all slopes (uslp, vslp,
C -------------------------------		 wslpi,wslpj )
#ifdef key_mpp
C
C ... Mpp: export boundary values to neighbouring processors
CCC 25/06/2003 Cristiano - It has been adjusted (add a parameter), but I am not sure of the value
CCC      CALL mpplnk( uslp , 2  )
CCC      CALL mpplnk( vslp , 3  )
CCC      CALL mpplnk( wslpi, 1  )
CCC      CALL mpplnk( wslpj, 1  )
      CALL mpplnk( uslp , 2, 1)
      CALL mpplnk( vslp , 3, 1)
      CALL mpplnk( wslpi, 1, 1)
      CALL mpplnk( wslpj, 1, 1)

#  else
C
C ... mono- or macro-tasking : U-, V- and T-points, 3D arrays, jk-slab
      CALL lbc( uslp , 2, 0, 1, ktask+1, jpkm1, 1 )
      CALL lbc( vslp , 3, 0, 1, ktask+1, jpkm1, 1 )
      CALL lbc( wslpi, 1, 0, 1, ktask+1, jpkm1, 1 )
      CALL lbc( wslpj, 1, 0, 1, ktask+1, jpkm1, 1 )
#endif

#  else

C No isopycnal diffusion used

#endif

      RETURN
      END
