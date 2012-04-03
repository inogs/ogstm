       MODULE stafun
CCC---------------------------------------------------------------------
CCC
CCC                         STAFUN
CCC                       **********
CCC
CCC  PURPOSE :
CCC  ---------
CCC     Statement function file: to be include in all routines




!      INTEGER kfi, kfj, kfk
!      REAL(8) pfs, pfp
!      REAL(8) pfb, pfn, pfa
!      REAL(8) pfx1, pfx2, pfy1, pfy2, pfz1, pfz2, pfu, pfv, pfw
!      REAL(8) fsfzpt
!      REAL(8) fsass
!      REAL(8) fsx, fsy, fsz
!      REAL(8) fsdept, fsdepu, fsdepv, fsdepf, fsdepw, fsdepuw, fsdepvw
!      REAL(8) fse3t, fse3u, fse3v, fse3f, fse3w, fse3uw, fse3vw
!      REAL(8) fsahtu, fsahtv, fsahtw, fsahtt
!#if defined key_trahdfeiv
!      REAL(8) fsaeiu, fsaeiv, fsaeiw
!#endif
!#if defined key_s_coord
!      REAL(8) fsde3w
!#endif
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

      USE myalloc
      IMPLICIT NONE

      PUBLIC
      CONTAINS


      REAL(8) FUNCTION fsfzpt( pfs, pfp )
      IMPLICIT NONE
      REAL(8), INTENT(IN) ::   pfs, pfp
      fsfzpt = ( -0.0575 + 1.710523e-3 * sqrt(pfs)  - 2.154996e-4 *pfs)*pfs - 7.53e-4 * pfp
      END FUNCTION fsfzpt
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
      REAL(8) FUNCTION fsass( pfb, pfn, pfa )
      IMPLICIT NONE
      REAL(8), INTENT(IN) ::  pfb, pfn, pfa
      fsass = atfp * ( pfb + pfa ) + atfp1 * pfn
      END FUNCTION
C
C-----------------------------------------------------------------------
C
C Up Stream Advection Scheme
C ==========================
C  fsx: along i-direction
C  fsy: along j-direction
C  fsz: along k-direction
C
      REAL(8) FUNCTION fsx( pfx1, pfx2, pfu )
      IMPLICIT NONE
      REAL(8), INTENT(IN) :: pfx1, pfx2, pfu
      fsx = ( ( pfu + abs(pfu) ) * pfx1+( pfu - abs(pfu) ) * pfx2 ) * 0.5
      END FUNCTION fsx

      REAL(8) FUNCTION fsy( pfy1, pfy2, pfv  )
      IMPLICIT NONE
      REAL(8), INTENT(IN) :: pfy1, pfy2, pfv
      fsy = ( ( pfv + abs(pfv) ) * pfy1 +( pfv - abs(pfv) ) * pfy2 ) * 0.5
      END FUNCTION fsy


      REAL(8) FUNCTION fsz( pfz1, pfz2, pfw )
      IMPLICIT NONE
      REAL(8), INTENT(IN) :: pfz1, pfz2, pfw
      fsz = ( ( pfw + abs(pfw) ) * pfz1+( pfw - abs(pfw) ) * pfz2 ) * 0.5
      END FUNCTION fsz


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
      REAL(8) FUNCTION fsdept ( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsdept = hbatt(kfi,kfj) * gsigt(kfk)
      END FUNCTION fsdept

      REAL(8) FUNCTION fsdepu ( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsdepu = hbatu(kfi,kfj) * gsigt(kfk)
      END FUNCTION fsdepu

      REAL(8) FUNCTION fsdepv ( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsdepv  = hbatv(kfi,kfj) * gsigt(kfk)
      END FUNCTION fsdepv

      REAL(8) FUNCTION fsdepv ( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsdepf = hbatf(kfi,kfj) * gsigt(kfk)
      END FUNCTION fsdepf



      REAL(8) FUNCTION fsdepw ( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsdepw = hbatt(kfi,kfj) * gsigw(kfk)
      END FUNCTION fsdepw

      REAL(8) FUNCTION fsdepuw ( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsdepuw( kfi, kfj, kfk ) = hbatu(kfi,kfj) * gsi3w(kfk)
      END FUNCTION fsdepuw

      REAL(8) FUNCTION fsdepvw ( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsdepvw = hbatv(kfi,kfj) * gsi3w(kfk)
      END FUNCTION fsdepvw

      REAL(8) FUNCTION fsde3w ( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsde3w = hbatt(kfi,kfj) * gsi3w(kfk)
      END FUNCTION fsde3w

!     ******************************************************

      REAL(8) FUNCTION fse3u( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fse3u = hbatu(kfi,kfj) * esigt(kfk)
      END FUNCTION fse3u

      REAL(8) FUNCTION fse3v( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fse3v  = hbatv(kfi,kfj) * esigt(kfk)
      END FUNCTION fse3v

      REAL(8) FUNCTION fse3f( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fse3f  = hbatf(kfi,kfj) * esigt(kfk)
      END FUNCTION fse3f

      REAL(8) FUNCTION fse3w( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fse3w ( kfi, kfj, kfk ) = hbatt(kfi,kfj) * esigw(kfk)
      END FUNCTION fse3w

      REAL(8) FUNCTION fse3uw( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fse3uw = hbatu(kfi,kfj) * esigw(kfk)
      END FUNCTION fse3uw

      REAL(8) FUNCTION fse3vw( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fse3vw = hbatv(kfi,kfj) * esigw(kfk)
      END FUNCTION fse3vw


C
# else

      REAL(8) FUNCTION fsdept( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsdept = gdept(kfk)
      END FUNCTION fsdept

      REAL(8) FUNCTION fsdepw( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsdepw = gdepw(kfk)
      END FUNCTION fsdepw

      REAL(8) FUNCTION fse3t( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fse3t = e3t(kfk)
      END FUNCTION fse3t

      REAL(8) FUNCTION fse3u( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fse3u = e3t(kfk)
      END FUNCTION fse3u

      REAL(8) FUNCTION fse3v( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fse3v = e3t(kfk)
      END FUNCTION fse3v

      REAL(8) FUNCTION fse3f( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fse3f = e3t(kfk)
      END FUNCTION fse3f

      REAL(8) FUNCTION fse3w( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fse3w = e3w(kfk)
      END FUNCTION fse3w

      REAL(8) FUNCTION fse3uw( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fse3uw = e3w(kfk)
      END FUNCTION fse3uw

      REAL(8) FUNCTION fse3vw( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fse3vw = e3w(kfk)
      END FUNCTION fse3vw



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
      REAL(8) FUNCTION fsahtt( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsahtt = ahtt(kfi,kfj,kfk)
      END FUNCTION fsahtt
#   else
      REAL(8) FUNCTION fsahtu( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsahtu = ahtu(kfi,kfj,kfk)
      END FUNCTION  fsahtu

      REAL(8) FUNCTION fsahtv( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsahtv = ahtv(kfi,kfj,kfk)
      END FUNCTION fsahtv

#    if defined key_trahdfiso || defined key_trahdfgeop
      REAL(8) FUNCTION fsahtw( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsahtw = ahtw(kfi,kfj,kfk)
      END FUNCTION fsahtw
#    endif
#  endif
# elif defined key_trahdfcoef2d
C 2D coefficient
#  ifdef key_trahdfbilap
      REAL(8) FUNCTION fsahtt( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsahtt = ahtt(kfi,kfj)
      END FUNCTION fsahtt
#   else
      REAL(8) FUNCTION fsahtu( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsahtu = ahtu(kfi,kfj)
      END FUNCTION fsahtu

      REAL(8) FUNCTION fsahtv( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsahtv = ahtv(kfi,kfj)
      END FUNCTION fsahtv

      REAL(8) FUNCTION fsahtw( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsahtw= ahtw(kfi,kfj)
      END FUNCTION fsahtw
#  endif
# elif defined key_trahdfcoef1d
C 1D coefficient
#  ifdef key_trahdfbilap
      REAL(8) FUNCTION fsahtt( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsahtt = ahtt(kfk)
      END FUNCTION fsahtt
#   else
      REAL(8) FUNCTION fsahtu( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsahtu = ahtu(kfk)
      END FUNCTION fsahtu

      REAL(8) FUNCTION fsahtv( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsahtv = ahtv(kfk)
      END FUNCTION fsahtv

      REAL(8) FUNCTION fsahtw( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsahtw = ahtw(kfk)
      END FUNCTION fsahtw
#  endif
# else
#  if defined key_off_degrad
C Constant coefficient but because of degradation...
      REAL(8) FUNCTION fsahtt( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsahtt= aht0
      END FUNCTION fsahtt

      REAL(8) FUNCTION fsahtu( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsahtu = ahtu( kfi, kfj, kfk )
      END FUNCTION fsahtu

      REAL(8) FUNCTION fsahtv( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsahtv = ahtv( kfi, kfj, kfk )
      END FUNCTION fsahtv

      REAL(8) FUNCTION fsahtw( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsahtw = ahtw( kfi, kfj, kfk )
      END FUNCTION fsahtw


#  else
C Constant coefficient
      REAL(8) FUNCTION fsahtt( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsahtt = aht0
      END FUNCTION fsahtt

      REAL(8) FUNCTION fsahtu( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsahtu = aht0
      END FUNCTION fsahtu

      REAL(8) FUNCTION fsahtv( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsahtv = aht0
      END FUNCTION fsahtv

      REAL(8) FUNCTION fsahtw( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsahtw = aht0
      END FUNCTION fsahtw

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
      REAL(8) FUNCTION fsaeiu( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsaeiu = aeiu(kfi,kfj,kfk)
      END FUNCTION fsaeiu

      REAL(8) FUNCTION fsaeiv( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsaeiv = aeiv(kfi,kfj,kfk)
      END FUNCTION fsaeiv

      REAL(8) FUNCTION fsaeiw( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsaeiw = aeiw(kfi,kfj,kfk)
      END FUNCTION fsaeiw


#  elif defined key_trahdfcoef2d
C 2D coefficient
      REAL(8) FUNCTION fsaeiu( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsaeiu = aeiu(kfi,kfj)
      END FUNCTION fsaeiu

      REAL(8) FUNCTION fsaeiv( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsaeiv = aeiv(kfi,kfj)
      END FUNCTION fsaeiv

      REAL(8) FUNCTION fsaeiw( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsaeiw = aeiw(kfi,kfj)
      END FUNCTION fsaeiw


#  elif defined key_trahdfcoef1d
C 1D coefficient
      REAL(8) FUNCTION fsaeiu( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsaeiu = aeiu(kfk)
      END FUNCTION fsaeiu

      REAL(8) FUNCTION fsaeiv( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsaeiv = aeiv(kfk)
      END FUNCTION fsaeiv

      REAL(8) FUNCTION fsaeiw( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsaeiw = aeiw(kfk)
      END FUNCTION fsaeiw

#  else
#  if defined key_off_degrad
C Constant coefficient but because of degradation...
      REAL(8) FUNCTION fsaeiu( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsaeiu = aeiu( kfi, kfj, kfk )
      END FUNCTION fsaeiu

      REAL(8) FUNCTION fsaeiv( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsaeiv = aeiv( kfi, kfj, kfk )
      END FUNCTION fsaeiv

      REAL(8) FUNCTION fsaeiw( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsaeiw = aeiw( kfi, kfj, kfk )
      END FUNCTION fsaeiw

#  else
C Constant coefficient
      REAL(8) FUNCTION fsaeiu( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsaeiu = aeiv0
      END FUNCTION fsaeiu

      REAL(8) FUNCTION fsaeiv( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsaeiv= aeiv0
      END FUNCTION fsaeiv

      REAL(8) FUNCTION fsaeiw( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsaeiw = aeiv0
      END FUNCTION fsaeiw


# endif
# endif
C
#endif
C
CC-----------------------------------------------------------------------

      END MODULE stafun
