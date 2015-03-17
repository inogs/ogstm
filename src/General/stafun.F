       MODULE stafun
CCC---------------------------------------------------------------------
CCC
CCC                         STAFUN
CCC                       **********
CCC
CCC  PURPOSE :
CCC  ---------
CCC     Statement function file: to be include in all routines




C-----------------------------------------------------------------------

C
C Up Stream Advection Scheme
C ==========================
C  fsx: along i-direction
C  fsy: along j-direction
C  fsz: along k-direction
      USE myalloc
      IMPLICIT NONE

      PUBLIC

      CONTAINS
C
      REAL(8) FUNCTION fsx( pfx1, pfx2, pfu )
      IMPLICIT NONE
      REAL(8), INTENT(IN) :: pfx1, pfx2, pfu
      REAL(8) abspfu
      abspfu = abs(pfu)
      fsx = ( ( pfu + abspfu ) * pfx1+( pfu - abspfu ) * pfx2 ) * 0.5
      END FUNCTION fsx

      REAL(8) FUNCTION fsy( pfy1, pfy2, pfv  )
      IMPLICIT NONE
      REAL(8), INTENT(IN) :: pfy1, pfy2, pfv
      REAL(8) abspfv
      abspfv = abs(pfv)
      fsy = ( ( pfv + abspfv ) * pfy1 +( pfv - abspfv ) * pfy2 ) * 0.5
      END FUNCTION fsy


      REAL(8) FUNCTION fsz( pfz1, pfz2, pfw )
      IMPLICIT NONE
      REAL(8), INTENT(IN) :: pfz1, pfz2, pfw
      REAL(8) abspfw
      abspfw = abs(pfw)
      fsz = ( ( pfw + abspfw ) * pfz1+( pfw - abspfw ) * pfz2 ) * 0.5
      END FUNCTION fsz


C
C-----------------------------------------------------------------------
C
C Vertical mesh
C =============
C  z-coordinates (default option) depth and vertical scale factors are
C                defined from 1d fields;
C  fsdept, fsdepw  : depth of model level at t- and w-points
C  fse3t, fse3u,   : factors at t-, u-, v-, f-, w-,  uw-, vw-points
C  fse3v, fse3f,     z-coordinates (default option) 
C  fse3w, fse3uw,
C  fse3vw 
C  (kfi,kfj,kfk)   : indexes of the position
C
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

C
CC-----------------------------------------------------------------------

      END MODULE stafun
