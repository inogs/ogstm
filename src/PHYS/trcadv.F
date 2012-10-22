C $Id: trcadv.F,v 1.2 2009-09-11 09:20:56 cvsogs01 Exp $

      SUBROUTINE trcadv(ktask,kt)
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE trcadv
CCC                     ******************
CCC
CCC  Purpose :
CCC  ---------
CCC     Compute the now trend due to the advection of tracers
CCC     (tr) and add it to the general trend of passive tracer equations.
CCC
CCC     Differents schemes are possible depending on the key
CCC
CCC     'key_trc_arakawa'
CCC     'key_trc_smolar'
CCC     'key_trc_tvd'
CCC     'key_trc_muscl'
CCC     'key_trc_ppm'
CCC

CC   MODIFICATIONS:
CC   --------------
CC      original  : 00-10 M.A. Foujols and E. Kestenare
CC----------------------------------------------------------------------
CC parameters and commons
CC ======================

       USE myalloc
       USE myalloc_mpp
       USE ADV_mem
       USE stafun

        IMPLICIT NONE

CC local declarations
CC ==================
      INTEGER ktask,kt


#if defined key_passivetrc
#    if defined key_trc_arakawa
#       include "trcadv.arakawa.h"
#    elif defined key_trc_smolar
         call smolar
!#       include "trcadv.smolar.h"
#    elif defined key_trc_upstream
#       include "trcadv.upstream.h"
#    elif defined key_trc_tvd
#       include "trcadv.tvd.h"
#    elif defined key_trc_muscl
#       include "trcadv.muscl.h"
#    elif defined key_trc_ppm
#       include "trcadv.ppm.h"
#    endif
C
C
#else
C
C no passive tracers
C
#endif
C
C
      RETURN
      END
