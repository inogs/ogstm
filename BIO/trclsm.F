C $Id: trclsm.F,v 1.2 2009-09-11 09:20:56 cvsogs01 Exp $

      SUBROUTINE trclsm
CCC---------------------------------------------------------------------


CCC
CCC                       ROUTINE trclsm
CCC                     ******************
CCC
CCC  PURPOSE :
CCC  ---------
CCC     READS specific NAMELIST for sms terms 
CCC
CC   METHOD :                   : no
CC   -------
CC
CC   INPUT :
CC   -----
CC
CC   OUTPUT :
CC   ------
CC
CC   WORKSPACE :                : no
CC   ---------
CC
CC   MODIFICATIONS:
CC   --------------
CC      original  : 99-10 (M.A. Foujols, M. Levy) passive tracer
CC                  04/00 (O. Aumont, M.A. Foujols) HAMOCC3 and P3ZD
CC----------------------------------------------------------------------
CC parameters and commons
CC ======================
c#include "parameter.h"
c#include "common.h"
       USE myalloc
CC+CC  Implicit typing is never allowed
        IMPLICIT NONE
CC+CC  Implicit typing is never allowed


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
C
C no passive tracers
C
#endif
      RETURN
      END
