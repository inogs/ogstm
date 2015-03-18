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


C
CC----------------------------------------------------------------------
C
C Lateral eddy diffusivity coefficient for passive tracer:
C ========================================================
C

      CONTAINS

C 1D coefficient
#  ifdef key_trc_hdfbilap
#   else
#  endif


      END MODULE STAFUN_PASSTRC
