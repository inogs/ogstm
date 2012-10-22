      MODULE STAFUN_PASSTRC
       USE myalloc
CCC
CCC                         STAFUN.TRC
CCC                       *************
CCC
CCC  PURPOSE :
CCC  ---------
CCC     Statement function file: to be include in all routines
CCC                  concerning passive tracer model


!      REAL(8) fsahtru, fsahtrv, fsahtrw, fsahtrt
!#if defined key_trc_hdfeiv
!      REAL(8) fsaeitru, fsaeitrv, fsaeitrw
!#endif
C
CC----------------------------------------------------------------------
C
C Lateral eddy diffusivity coefficient for passive tracer:
C ========================================================
C  fsahtru, fsahtrv, : lateral eddy diffusivity coef. at u-, v-, w-points
C  fsahtrw             (for second order diffusive operator on tracers)
C  fsahtrt           : lateral eddy diffusivity coef. at t-point
C                    (for fourth order diffusive operator on tracers)
C  (kfi,kfj,kfk)     :  indexes of the position
C

      CONTAINS
#ifdef key_trahdfcoef3d
C 3D coefficient
#  ifdef key_trc_hdfbilap
      REAL(8) FUNCTION fsahtrt( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsahtrt = trcrat * ahtt(kfi,kfj,kfk)
      END FUNCTION fsahtrt
#   else
      REAL(8) FUNCTION fsahtru( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsahtru = trcrat * ahtu(kfi,kfj,kfk)
      END FUNCTION fsahtru

      REAL(8) FUNCTION fsahtrv( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsahtrv = trcrat * ahtv(kfi,kfj,kfk)
      END FUNCTION fsahtrv


#    if defined key_trc_hdfiso || defined key_trc_hdfgeop
      REAL(8) FUNCTION fsahtrw( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsahtrw = trcrat * ahtw(kfi,kfj,kfk)
      END FUNCTION fsahtrw
#    endif
#  endif
# elif defined key_trahdfcoef2d
C 2D coefficient
#  ifdef key_trc_hdfbilap
      REAL(8) FUNCTION fsahtrt( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsahtrt = trcrat * ahtt(kfi,kfj)
      END FUNCTION fsahtrt
#   else
      REAL(8) FUNCTION fsahtru( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsahtru =  trcrat * ahtu(kfi,kfj)
      END FUNCTION fsahtru

      REAL(8) FUNCTION fsahtrv( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsahtrv =  trcrat * ahtv(kfi,kfj)
      FUNCTION fsahtrv

      REAL(8) FUNCTION fsahtrw( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsahtrw =  trcrat * ahtw(kfi,kfj)
      FUNCTION fsahtrw
#  endif
# elif defined key_trahdfcoef1d
C 1D coefficient
#  ifdef key_trc_hdfbilap
      REAL(8) FUNCTION fsahtrt( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsahtrt = trcrat * ahtt(kfk)
      END FUNCTION fsahtrt
#   else
      REAL(8) FUNCTION fsahtru( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsahtru = trcrat * ahtu(kfk)
      END FUNCTION fsahtru

      REAL(8) FUNCTION fsahtrv( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsahtrv = trcrat * ahtv(kfk)
      END FUNCTION fsahtrv

      REAL(8) FUNCTION fsahtrw( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsahtrw = trcrat * ahtw(kfk)
      END FUNCTION fsahtrw
#  endif
# else
C Constant coefficient
      REAL(8) FUNCTION fsahtrt( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsahtrt = ahtrc0
      END FUNCTION fsahtrt

      REAL(8) FUNCTION fsahtru( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsahtru = ahtrc0
      END FUNCTION fsahtru

      REAL(8) FUNCTION fsahtrv( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsahtrv = ahtrc0
      END FUNCTION fsahtrv

      REAL(8) FUNCTION fsahtrw( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsahtrw = ahtrc0
      END FUNCTION fsahtrw


#endif
C
#if defined key_trc_hdfeiv
CC----------------------------------------------------------------------
C
C Eddy induced velocity coefficient for passive tracer:
C ====================================================
C  fsaeiu, fsaeiv, : eddy induced velocity coefficients at u-, v- and
C  fsaeiw            w-points 
C  (kfi,kfj,kfk)   :  indexes of the position
C
# ifdef key_trahdfcoef3d
C 3D coefficient
      REAL(8) FUNCTION fsaeitru( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsaeitru =  trcrat * aeiu(kfi,kfj,kfk)
      END FUNCTION fsaeitru

      REAL(8) FUNCTION fsaeitrv( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsaeitrv=  trcrat * aeiv(kfi,kfj,kfk)
      END FUNCTION fsaeitrv

      REAL(8) FUNCTION fsaeitrw( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsaeitrw =  trcrat * aeiw(kfi,kfj,kfk)
      END FUNCTION fsaeitrw


#  elif defined key_trahdfcoef2d
C 2D coefficient

      REAL(8) FUNCTION fsaeitru( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsaeitru =  trcrat * aeiu(kfi,kfj)
      END FUNCTION fsaeitru

      REAL(8) FUNCTION fsaeitrv( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsaeitrv =  trcrat * aeiv(kfi,kfj)
      END FUNCTION fsaeitrv

      REAL(8) FUNCTION fsaeitrw( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsaeitrw =  trcrat * aeiw(kfi,kfj)
      END FUNCTION fsaeitrw


#  elif defined key_trahdfcoef1d
C 1D coefficient

      REAL(8) FUNCTION fsaeitru( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsaeitru =  trcrat * aeiu(kfk)
      END FUNCTION fsaeitru

      REAL(8) FUNCTION fsaeitrv( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsaeitrv =  trcrat * aeiv(kfk)
      END FUNCTION fsaeitrv

      REAL(8) FUNCTION fsaeitrw( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsaeitrw =  trcrat * aeiw(kfk)
      END FUNCTION fsaeitrw



#  else
C Constant coefficient
      REAL(8) FUNCTION fsaeitru( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsaeitru = aeivtr0
      END FUNCTION fsaeitru

      REAL(8) FUNCTION fsaeitrv( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsaeitrv = aeivtr0
      END FUNCTION fsaeitrv

      REAL(8) FUNCTION fsaeitrw( kfi, kfj, kfk )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  kfi, kfj, kfk
      fsaeitrw = aeivtr0
      END FUNCTION fsaeitrw


# endif

#endif

      END MODULE STAFUN_PASSTRC
