CCC $Header: /cvsroot/opatm-bfm/opa_model/OPA/stafun.passivetrc.h,v 1.2 2009-09-11 09:20:56 cvsogs01 Exp $
CCC---------------------------------------------------------------------
CCC
CCC                         STAFUN.TRC
CCC                       *************
CCC
CCC  PURPOSE :
CCC  ---------
CCC     Statement function file: to be include in all routines
CCC                  concerning passive tracer model
CCC
CC   MODIFICATIONS :
CC   -------------
CC      original : 00-11 (M.A. Foujols E. Kestenare)
CC      addition :
CCC---------------------------------------------------------------------
CCC  OPA8, LODYC (2000)
CCC---------------------------------------------------------------------
      REAL(8) fsahtru, fsahtrv, fsahtrw, fsahtrt
#if defined key_trc_hdfeiv
      REAL(8) fsaeitru, fsaeitrv, fsaeitrw
#endif
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
#ifdef key_trahdfcoef3d
C 3D coefficient
#  ifdef key_trc_hdfbilap
      fsahtrt( kfi, kfj, kfk ) = trcrat * ahtt(kfi,kfj,kfk)
#   else
      fsahtru( kfi, kfj, kfk ) = trcrat * ahtu(kfi,kfj,kfk)
      fsahtrv( kfi, kfj, kfk ) = trcrat * ahtv(kfi,kfj,kfk)
#    if defined key_trc_hdfiso || defined key_trc_hdfgeop
      fsahtrw( kfi, kfj, kfk ) = trcrat * ahtw(kfi,kfj,kfk)
#    endif
#  endif
# elif defined key_trahdfcoef2d
C 2D coefficient
#  ifdef key_trc_hdfbilap
      fsahtrt( kfi, kfj, kfk ) = trcrat * ahtt(kfi,kfj)
#   else
      fsahtru( kfi, kfj, kfk ) =  trcrat * ahtu(kfi,kfj)
      fsahtrv( kfi, kfj, kfk ) =  trcrat * ahtv(kfi,kfj)
      fsahtrw( kfi, kfj, kfk ) =  trcrat * ahtw(kfi,kfj)
#  endif
# elif defined key_trahdfcoef1d
C 1D coefficient
#  ifdef key_trc_hdfbilap
      fsahtrt( kfi, kfj, kfk ) = trcrat * ahtt(kfk)
#   else
      fsahtru( kfi, kfj, kfk ) = trcrat * ahtu(kfk)
      fsahtrv( kfi, kfj, kfk ) = trcrat * ahtv(kfk)
      fsahtrw( kfi, kfj, kfk ) = trcrat * ahtw(kfk)
#  endif
# else
C Constant coefficient
      fsahtrt( kfi, kfj, kfk ) = ahtrc0
      fsahtru( kfi, kfj, kfk ) = ahtrc0
      fsahtrv( kfi, kfj, kfk ) = ahtrc0
      fsahtrw( kfi, kfj, kfk ) = ahtrc0
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
      fsaeitru( kfi, kfj, kfk ) =  trcrat * aeiu(kfi,kfj,kfk)
      fsaeitrv( kfi, kfj, kfk ) =  trcrat * aeiv(kfi,kfj,kfk)
      fsaeitrw( kfi, kfj, kfk ) =  trcrat * aeiw(kfi,kfj,kfk)
#  elif defined key_trahdfcoef2d
C 2D coefficient
      fsaeitru( kfi, kfj, kfk ) =  trcrat * aeiu(kfi,kfj)
      fsaeitrv( kfi, kfj, kfk ) =  trcrat * aeiv(kfi,kfj)
      fsaeitrw( kfi, kfj, kfk ) =  trcrat * aeiw(kfi,kfj)
#  elif defined key_trahdfcoef1d
C 1D coefficient
      fsaeitru( kfi, kfj, kfk ) =  trcrat * aeiu(kfk)
      fsaeitrv( kfi, kfj, kfk ) =  trcrat * aeiv(kfk)
      fsaeitrw( kfi, kfj, kfk ) =  trcrat * aeiw(kfk)
#  else
C Constant coefficient
      fsaeitru( kfi, kfj, kfk ) = aeivtr0
      fsaeitrv( kfi, kfj, kfk ) = aeivtr0
      fsaeitrw( kfi, kfj, kfk ) = aeivtr0
# endif
C
#endif
C
CC----------------------------------------------------------------------
