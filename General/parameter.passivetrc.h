c $Id: parameter.passivetrc.h,v 1.2 2009-09-11 09:20:56 cvsogs01 Exp $
#if defined key_passivetrc
CCC---------------------------------------------------------------------
CCC
CCC                         PARAMETER passivetrc
CCC                       ************************
CCC
CCC  PURPOSE :
CCC  ---------
CCC     Include parameter FILE for passive tracer
CCC
CC   MODIFICATIONS :
CC   -------------
CC	original : 96 (M. Levy)
CC                 07/99 (M. Levy for NNPZDDOM or NPZD model)
CC                 04/00 (O. Aumont, M.A. Foujols) HAMOCC3 and P3ZD
CC
CCC---------------------------------------------------------------------
CCC  OPA8, LODYC (15/11/96)
CCC---------------------------------------------------------------------
CC
CC passive tracers
CC ---------------
CC       jptra  : number of passive tracers
CC
#    if defined key_trc_npzd
#    include "parameter.passivetrc.npzd.h"
#    elif defined key_trc_bfm
#    include "parameter.passivetrc.bfm.h"
#    elif defined key_trc_npzdb
#    include "parameter.passivetrc.npzdb.h"
#    elif defined key_trc_generic
#    include "parameter.passivetrc.generic.h"
#    else
CC    default CASE : temperature and salinity as passive tracers
      INTEGER jptra
      PARAMETER (jptra = 2)
#    endif
#    if defined key_trc_diatrd
CC
CC number of dynamical trends
CC --------------------------
CC
      INTEGER jpdiatrc
#        if defined key_trahdfeiv
CC
CC we keep 3 more trends for eddy induced flux
CC (gent velocity)
CC
#            if defined key_trc_dmp
      PARAMETER (jpdiatrc = 10)
#            else
      PARAMETER (jpdiatrc = 9)
#            endif
#        else
#            if defined key_trc_dmp
      PARAMETER (jpdiatrc = 7)
#            else
      PARAMETER (jpdiatrc = 6)
#            endif
#        endif
#    endif
#    if defined key_trc_diaadd
CC
CC possibility for additional 3d and 2d output
CC -------------------------------------------
CC
      INTEGER jpdia3d, jpdia2d
      PARAMETER (jpdia3d = 1, jpdia2d = 13)
#    endif
#else
CC
CC no passive tracer 
CC
#endif
