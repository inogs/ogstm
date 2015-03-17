      SUBROUTINE inihdf
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE inihdf
CCC                     ******************
CCC
CCC  Purpose :
CCC  ---------
CCC     Initializations of the horizontal coefficients
CCC
CC   Method :
CC   -------
CC
CC    Eddy diffusivity and eddy induced velocity cefficients:
CC         default option   : constant coef. aht0 (namelist)
CC        'key_trahdfcoef1d': depth dependent coef. defined in 
CC                            in 'inihdf.tra.coef3d.h' file
CC
CC    Direction of lateral diffusion (tracers and/or momentum)
CC
CC    Biharmonic operator:
CC
CC   ---------
CC
CC   Reference :
CC   ----------
CC      Madec, G. and M. Imbard, 1996, A global ocean mesh to overcome
CC      the North Pole singularity, Climate Dynamics, 12, 381-388.
       USE myalloc
       USE stafun
        IMPLICIT NONE
CC ==================
      INTEGER  jk
C

      REAL(8) zkah, zahr, za00, za01, zahf, zahs
      REAL(8) zahtf, zahts


C 1. Lateral eddy diffusivity and eddy induced velocity coefficients
C ==================================================================
C     user: modify the 'inihdf.tra.coef_d.h' include files following
C           your desiderata

C
C 1D coefficients : 
C ----------------
#     include "inihdf.tra.coef1d.h"
C


      END SUBROUTINE inihdf
