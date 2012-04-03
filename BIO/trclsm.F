
      SUBROUTINE trclsm
CCC---------------------------------------------------------------------


CCC
CCC                       ROUTINE trclsm
CCC                     ******************
CCC
CCC  PURPOSE :
CCC  ---------
CCC     READS specific NAMELIST for sms terms 
CC
CC   MODIFICATIONS:
CC   --------------
CC      original  : 99-10 (M.A. Foujols, M. Levy) passive tracer
CC                  04/00 (O. Aumont, M.A. Foujols) HAMOCC3 and P3ZD
CC----------------------------------------------------------------------
CC parameters and commons
CC ======================

       USE myalloc
       USE stafun
        IMPLICIT NONE



#if defined key_passivetrc
#    if defined key_trc_npzd
#    include "trclsm.npzd.h"
#    elif defined key_trc_npzdb
#    include "trclsm.opt.h"
#    include "trclsm.sed.h"
#    elif defined key_trc_bfm
#    include "trclsm.opt.h"
#    include "trclsm.sed.h"
#    elif defined key_trc_generic
#    else
C      user can define here his own sms namelist
#    endif
#else

C no passive tracers

#endif

      END SUBROUTINE trclsm
