CC $Header: /cvsroot/opatm-bfm/opa_model/OPA/stafun.h,v 1.2 2009-09-11 09:20:56 cvsogs01 Exp $
CCC---------------------------------------------------------------------
CCC
CCC                         STAFUN
CCC                       **********
CCC
CCC  PURPOSE :
CCC  ---------
CCC     Statement function file: to be include in all routines
CCC
CC   MODIFICATIONS :
CC   -------------
CC      original : 95-12 (G. Madec)
CC      addition :
CCC---------------------------------------------------------------------
CCC  OPA8, LODYC (1997)
CCC---------------------------------------------------------------------
      INTEGER kfi, kfj, kfk
      REAL(8) pfs, pfp
      REAL(8) pfb, pfn, pfa
      REAL(8) pfx1, pfx2, pfy1, pfy2, pfz1, pfz2, pfu, pfv, pfw
      REAL(8) fsfzpt
      REAL(8) fsass
      REAL(8) fsx, fsy, fsz
      REAL(8) fsdept, fsdepu, fsdepv, fsdepf, fsdepw, fsdepuw, fsdepvw
      REAL(8) fse3t, fse3u, fse3v, fse3f, fse3w, fse3uw, fse3vw
      REAL(8) fsahtu, fsahtv, fsahtw, fsahtt
#if defined key_trahdfeiv
      REAL(8) fsaeiu, fsaeiv, fsaeiw
#endif
#if defined key_s_coord
      REAL(8) fsde3w
#endif
C-----------------------------------------------------------------------
C
C Ice freezing point
C ==================
C  fsfzpt: freezing point of seawater in degrees celsius
C       units : salinity        pfs       (ipss-78)
C               pressure        pfp      decibars
C               temperature     fszfpt   degrees celsius
C               freezing pt
C	reference : unesco tech. papers in the marine science no 28 1978
C               eigth report jpots
C               annex 6 freezing point of seawater F.J.Millero pp.29-35
C	checkvalue: fsfzpt=-2.588567 deg.c,for s=40.0,p=500 decibars
C
      fsfzpt( pfs, pfp ) = ( -0.0575 + 1.710523e-3 * sqrt(pfs)
     $                               - 2.154996e-4 *      pfs  ) * pfs
     $                   - 7.53e-4 * pfp
C
C-----------------------------------------------------------------------
C
C Asselin time filter
C ===================
C  fsass: asselin time filter
C           pfb : previous variable (before)
C           pfn : present  variable (now)
C           pfa : next     variable (after)
C
      fsass( pfb, pfn, pfa ) = atfp * ( pfb + pfa ) + atfp1 * pfn
C
C-----------------------------------------------------------------------
C
C Up Stream Advection Scheme
C ==========================
C  fsx: along i-direction
C  fsy: along j-direction
C  fsz: along k-direction
C
      fsx( pfx1, pfx2, pfu ) = ( ( pfu + abs(pfu) ) * pfx1
     $                          +( pfu - abs(pfu) ) * pfx2 ) * 0.5
      fsy( pfy1, pfy2, pfv ) = ( ( pfv + abs(pfv) ) * pfy1
     $                          +( pfv - abs(pfv) ) * pfy2 ) * 0.5
      fsz( pfz1, pfz2, pfw ) = ( ( pfw + abs(pfw) ) * pfz1
     $                          +( pfw - abs(pfw) ) * pfz2 ) * 0.5
C
C-----------------------------------------------------------------------
C
C Vertical mesh
C =============
C  z-coordinates (default option) depth and vertical scale factors are
C                defined from 1d fields;
C  s-coordinates (key_s_coord defined) depth and vertical scale factors
C		 are the product of a bathymetry field by a 1d coef.
C  fsdept, fsdepw  : depth of model level at t- and w-points
C  fse3t, fse3u,   : factors at t-, u-, v-, f-, w-,  uw-, vw-points
C  fse3v, fse3f,     z-coordinates (default option) 
C  fse3w, fse3uw,
C  fse3vw 
C  (kfi,kfj,kfk)   : indexes of the position
C
#ifdef key_s_coord
      fsdept ( kfi, kfj, kfk ) = hbatt(kfi,kfj) * gsigt(kfk)
      fsdepu ( kfi, kfj, kfk ) = hbatu(kfi,kfj) * gsigt(kfk)
      fsdepv ( kfi, kfj, kfk ) = hbatv(kfi,kfj) * gsigt(kfk)
      fsdepf ( kfi, kfj, kfk ) = hbatf(kfi,kfj) * gsigt(kfk)
      fsdepw ( kfi, kfj, kfk ) = hbatt(kfi,kfj) * gsigw(kfk)
      fsdepuw( kfi, kfj, kfk ) = hbatu(kfi,kfj) * gsi3w(kfk)
      fsdepvw( kfi, kfj, kfk ) = hbatv(kfi,kfj) * gsi3w(kfk)
      fsde3w ( kfi, kfj, kfk ) = hbatt(kfi,kfj) * gsi3w(kfk)
C
      fse3t ( kfi, kfj, kfk ) = hbatt(kfi,kfj) * esigt(kfk)
      fse3u ( kfi, kfj, kfk ) = hbatu(kfi,kfj) * esigt(kfk)
      fse3v ( kfi, kfj, kfk ) = hbatv(kfi,kfj) * esigt(kfk)
      fse3f ( kfi, kfj, kfk ) = hbatf(kfi,kfj) * esigt(kfk)
      fse3w ( kfi, kfj, kfk ) = hbatt(kfi,kfj) * esigw(kfk)
      fse3uw( kfi, kfj, kfk ) = hbatu(kfi,kfj) * esigw(kfk)
      fse3vw( kfi, kfj, kfk ) = hbatv(kfi,kfj) * esigw(kfk)
C
# else
      fsdept( kfi, kfj, kfk ) = gdept(kfk)
      fsdepw( kfi, kfj, kfk ) = gdepw(kfk)
      fse3t ( kfi, kfj, kfk ) = e3t(kfk)
      fse3u ( kfi, kfj, kfk ) = e3t(kfk)
      fse3v ( kfi, kfj, kfk ) = e3t(kfk)
      fse3f ( kfi, kfj, kfk ) = e3t(kfk)
      fse3w ( kfi, kfj, kfk ) = e3w(kfk)
      fse3uw( kfi, kfj, kfk ) = e3w(kfk)
      fse3vw( kfi, kfj, kfk ) = e3w(kfk)
C
#endif
C
CC----------------------------------------------------------------------
C
C Lateral eddy diffusivity coefficient:
C ====================================
C  fsahtu, fsahtv, : lateral eddy diffusivity coef. at u-, v-, w-points
C  fsahtw            (for second order diffusive operator on tracers)
C  fsahtt          : lateral eddy diffusivity coef. at t-point
C                    (for fourth order diffusive operator on tracers)
C  (kfi,kfj,kfk)   :  indexes of the position
C
#ifdef key_trahdfcoef3d
C 3D coefficient
#  ifdef key_trahdfbilap
      fsahtt( kfi, kfj, kfk ) = ahtt(kfi,kfj,kfk)
#   else
      fsahtu( kfi, kfj, kfk ) = ahtu(kfi,kfj,kfk)
      fsahtv( kfi, kfj, kfk ) = ahtv(kfi,kfj,kfk)
#    if defined key_trahdfiso || defined key_trahdfgeop
      fsahtw( kfi, kfj, kfk ) = ahtw(kfi,kfj,kfk)
#    endif
#  endif
# elif defined key_trahdfcoef2d
C 2D coefficient
#  ifdef key_trahdfbilap
      fsahtt( kfi, kfj, kfk ) = ahtt(kfi,kfj)
#   else
      fsahtu( kfi, kfj, kfk ) = ahtu(kfi,kfj)
      fsahtv( kfi, kfj, kfk ) = ahtv(kfi,kfj)
      fsahtw( kfi, kfj, kfk ) = ahtw(kfi,kfj)
#  endif
# elif defined key_trahdfcoef1d
C 1D coefficient
#  ifdef key_trahdfbilap
      fsahtt( kfi, kfj, kfk ) = ahtt(kfk)
#   else
      fsahtu( kfi, kfj, kfk ) = ahtu(kfk)
      fsahtv( kfi, kfj, kfk ) = ahtv(kfk)
      fsahtw( kfi, kfj, kfk ) = ahtw(kfk)
#  endif
# else
#  if defined key_off_degrad
C Constant coefficient but because of degradation...
      fsahtt( kfi, kfj, kfk ) = aht0
      fsahtu( kfi, kfj, kfk ) = ahtu( kfi, kfj, kfk )
      fsahtv( kfi, kfj, kfk ) = ahtv( kfi, kfj, kfk )
      fsahtw( kfi, kfj, kfk ) = ahtw( kfi, kfj, kfk )
#  else
C Constant coefficient
      fsahtt( kfi, kfj, kfk ) = aht0
      fsahtu( kfi, kfj, kfk ) = aht0
      fsahtv( kfi, kfj, kfk ) = aht0
      fsahtw( kfi, kfj, kfk ) = aht0
#   endif
#endif
C
#if defined key_trahdfeiv
CC----------------------------------------------------------------------
C
C Eddy induced velocity coefficient:
C =================================
C  fsaeiu, fsaeiv, : eddy induced velocity coefficients at u-, v- and
C  fsaeiw            w-points
C  (kfi,kfj,kfk)   :  indexes of the position
C
# ifdef key_trahdfcoef3d
C 3D coefficient
      fsaeiu( kfi, kfj, kfk ) = aeiu(kfi,kfj,kfk)
      fsaeiv( kfi, kfj, kfk ) = aeiv(kfi,kfj,kfk)
      fsaeiw( kfi, kfj, kfk ) = aeiw(kfi,kfj,kfk)
#  elif defined key_trahdfcoef2d
C 2D coefficient
      fsaeiu( kfi, kfj, kfk ) = aeiu(kfi,kfj)
      fsaeiv( kfi, kfj, kfk ) = aeiv(kfi,kfj)
      fsaeiw( kfi, kfj, kfk ) = aeiw(kfi,kfj)
#  elif defined key_trahdfcoef1d
C 1D coefficient
      fsaeiu( kfi, kfj, kfk ) = aeiu(kfk)
      fsaeiv( kfi, kfj, kfk ) = aeiv(kfk)
      fsaeiw( kfi, kfj, kfk ) = aeiw(kfk)
#  else
#  if defined key_off_degrad
C Constant coefficient but because of degradation...
      fsaeiu( kfi, kfj, kfk ) = aeiu( kfi, kfj, kfk )
      fsaeiv( kfi, kfj, kfk ) = aeiv( kfi, kfj, kfk )
      fsaeiw( kfi, kfj, kfk ) = aeiw( kfi, kfj, kfk )
#  else
C Constant coefficient
      fsaeiu( kfi, kfj, kfk ) = aeiv0
      fsaeiv( kfi, kfj, kfk ) = aeiv0
      fsaeiw( kfi, kfj, kfk ) = aeiv0
# endif
# endif
C
#endif
C
CC-----------------------------------------------------------------------

