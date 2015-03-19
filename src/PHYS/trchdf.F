      SUBROUTINE trchdf
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE trchdf
CCC                     ******************
CCC
CCC  PURPOSE :
CCC  ---------
CCC     compute the before horizontal passive tracer diffusive trend
CCC     add it to the general trend of tracer equation.
CCC
CCC     harmonic operator (default option)
CCC    -----------------
CCC       z-coordinates (default option):
CCC             default key : operator acting along model level surfaces
CCC                (which are geopotential surfaces).
CCC        'key_trc_hdfiso' : operator acting along neutral surfaces
CCC                (rotation of the diffusive tensor) plus an eddy
CCC            induced advection if 'key_trc_hdfeiv' defined.
CCC       s-coordinates ('key_s_coord'):
CCC        default key : operator acting along model level surfaces
CCC            (which are not geopotential surfaces).
CCC        'key_trc_hdfgeop' : operator acting along geopotential
CCC            surfaces (rotation of the diffusive tensor).
CCC        'key_trc_hdfiso' : operator acting along neutral surfaces
CCC                     (rotation of the diffusive tensor) plus an eddy
CCC                     induced advection if 'key_trc_hdfeiv' defined.
CCC
CCC    biharmonic operator ('key_trc_hdfbilap')
CCC    -------------------
CCC       z-coordinates (default option):
CCC             default key : operator acting along model level surfaces
CCC                     (which are geopotential surfaces).
CCC        s-coordinates ('key_s_coord'):
CCC             default key : operator acting along model level surfaces
CCC                     (which are not geopotential surfaces).
CCC             'key_trc_hdfgeop' : operator acting along geopotential
CCC                     surfaces (rotation of the diffusive tensor).
CCC
CCC
CCC     notes: 
CCC            IF you want no horizontal diffusion you have to switch
CCC            the boolean lhdf to .FALSE.
CCC
CCC     MODIFICATION : 00-11 (M.A.Foujols E. Kestenare)
CCC              differents keys for passive tracer
CCC


#ifdef key_trc_hdfbilap
#include "trchdf.bilaplacian.h"
# endif

#if defined key_trc_hdflap
#include "trchdf.laplacian.h"
# endif

      END SUBROUTINE trchdf
