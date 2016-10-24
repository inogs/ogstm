      SUBROUTINE trchdf
!!!---------------------------------------------------------------------
!!!
!!!                       ROUTINE trchdf
!!!                     ******************
!!!
!!!  PURPOSE :
!!!  ---------
!!!     compute the before horizontal passive tracer diffusive trend
!!!     add it to the general trend of tracer equation.
!!!
!!!     harmonic operator (default option)
!!!    -----------------
!!!       z-coordinates (default option):
!!!             default key : operator acting along model level surfaces
!!!                (which are geopotential surfaces).
!!!        'key_trc_hdfiso' : operator acting along neutral surfaces
!!!                (rotation of the diffusive tensor) plus an eddy
!!!            induced advection if 'key_trc_hdfeiv' defined.
!!!       s-coordinates ('key_s_coord'):
!!!        default key : operator acting along model level surfaces
!!!            (which are not geopotential surfaces).
!!!        'key_trc_hdfgeop' : operator acting along geopotential
!!!            surfaces (rotation of the diffusive tensor).
!!!        'key_trc_hdfiso' : operator acting along neutral surfaces
!!!                     (rotation of the diffusive tensor) plus an eddy
!!!                     induced advection if 'key_trc_hdfeiv' defined.
!!!
!!!    biharmonic operator ('key_trc_hdfbilap')
!!!    -------------------
!!!       z-coordinates (default option):
!!!             default key : operator acting along model level surfaces
!!!                     (which are geopotential surfaces).
!!!        s-coordinates ('key_s_coord'):
!!!             default key : operator acting along model level surfaces
!!!                     (which are not geopotential surfaces).
!!!             'key_trc_hdfgeop' : operator acting along geopotential
!!!                     surfaces (rotation of the diffusive tensor).
!!!
!!!
!!!     notes: 
!!!            IF you want no horizontal diffusion you have to switch
!!!            the boolean lhdf to .FALSE.
!!!
!!!     MODIFICATION : 00-11 (M.A.Foujols E. Kestenare)
!!!              differents keys for passive tracer
!!!

USE ogstm_mpi_module

#ifdef key_trc_hdfbilap
#include "trchdf.bilaplacian.h"
# endif

#if defined key_trc_hdflap
#include "trchdf.laplacian.h"
# endif

      END SUBROUTINE trchdf
