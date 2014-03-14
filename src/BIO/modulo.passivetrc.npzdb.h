CC
CC      vsed      : sedimentation speed (NAMELIST)
CC      remine()  : remineralisation trend
CC
CC      tmaxr     : maximum coefficient for passive tracer damping (NAMELIST)
CC      tminr     : minimum coefficient for passive tracer damping (NAMELIST)
CC      remdmp()  : damping coefficient of passive tracers (depth dependant)
CC

      REAL(8) vsed,tmaxr,tminr
      REAL(8), allocatable :: remine(:,:,:)
      REAL(8), allocatable :: remdmp(:,:)

CC

CC
CC----------------------------------------------------------------------
CC
CC COMMON/cotopt/ : optical parameters
CC -----------------------------------
CC
CC      xze       : euphotic layer depth
CC      xpar      : par (photosynthetic available radiation)
CC      xkr0      : water coefficient absorption in red (NAMELIST)
CC      xkg0      : water coefficient absorption in green (NAMELIST)
CC      xkrp      : pigment coefficient absorption in red (NAMELIST)
CC      xkgp      : pigment coefficient absorption in green (NAMELIST)
CC      xlr       : exposant for pigment absorption in red (NAMELIST)
CC      xlg       : exposant for pigment absorption in green (NAMELIST)
CC      rpig      : chla/chla+phea ratio (NAMELIST)
CC
      REAL(8) xkr0,xkg0,xkrp,xkgp,xlr,xlg,rpig

      REAL(8), allocatable :: xze(:,:)
      REAL(8), allocatable :: xpar(:,:,:)
CC----------------------------------------------------------------------
CC
