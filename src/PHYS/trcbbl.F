      SUBROUTINE trcbbl
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE trcbbl
CCC                     ******************
CCC
CCC  Purpose :
CCC  --------
CCC     Compute the before tracer (t & s) trend associated with the
CCC    bottom boundary layer and add it to the general trend of tracer
CCC    equations. The bottom boundary layer is supposed to be a purely
CCC     diffusive bottom boundary layer.
CCC
CC
CC
CC   Method :
CC   -------
CC    When the product grad( rho) * grad(h) < 0 (where grad is a
CC    along bottom slope gradient) an additional lateral 2nd order
CC    diffusion along the bottom slope is added to the general
CC    tracer trend, otherwise the additional trend is set to 0.
CC      Second order operator (laplacian type) with variable coefficient
CC      computed as follow for temperature (idem on s): 
CC         difft = 1/(e1t*e2t*e3t) { di-1[ ahbt e2u*e3u/e1u di[ztb] ]
CC                                 + dj-1[ ahbt e1v*e3v/e2v dj[ztb] ] }
CC    where ztb is a 2D array: the bottom ocean te;perature and ahtb
CC    is a time and space varying diffusive coefficient defined by:
CC       ahbt = zahbp    if grad(rho).grad(h) < 0
CC            = 0.       otherwise.
CC    Note that grad(.) is the along bottom slope gradient. grad(rho)
CC    is evaluated using the local density (i.e. referenced at the
CC    local depth). Typical value of ahbt is 2000 m2/s (equivalent to
CC    a downslope velocity of 20 cm/s if the condition for slope
CC    convection is satified)
CC
CC      Add this before trend to the general trend (ta,sa) of the 
CC    botton ocean tracer point:
CC              ta = ta + difft
CC
CC      'key_diatrends' defined, the trends are saved for diagnostics.
CC
CC      not macro-tasked
CC
CC   Input :
CC   ------
CC
CC   Output :
CC   -------
CC
CC   References :                 NO
CC   -----------
CC    Beckmann, A., and R. Doscher, 1997, J. Phys.Oceanogr., 581-591.

       USE myalloc
       USE stafun
        IMPLICIT NONE

#include <mpif.h>


#if defined key_passivetrc

#     if defined key_off_tra
#     include "trcbbl.offline.h"
#     else
#     include "trcbbl.online.h"
#     endif

#endif

      END SUBROUTINE trcbbl
