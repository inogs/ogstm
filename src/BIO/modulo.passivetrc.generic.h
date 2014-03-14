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

