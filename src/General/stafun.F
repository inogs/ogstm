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

      END MODULE stafun
