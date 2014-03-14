CCC---------------------------------------------------------------------
CCC
CCC                         PARAMETER passivetrc.bfm
CCC                       *****************************
CCC
CCC  purpose :
CCC  ---------
CCC     INCLUDE PARAMETER FILE for passive tracer generic model
CCC
CCC  --------------
CC
CC       jptra  : number of tracers
CC
      INTEGER jptra
      PARAMETER (jptra = 5)
CC
CC productive layer depth
CC ----------------------
CC       jpkb   : first vertical layers where generic tracer is active
CC       jpkbm1 : jpkb - 1
CC
      INTEGER jpkb,jpkbm1
      PARAMETER (jpkb = 32,jpkbm1 = 31)
CC
CC number of generic model trends
CC ---------------------------
CC
      INTEGER jpdiabio
      PARAMETER (jpdiabio = -1)
CC
CC    NOW ASSIGN A PARAMETER TO NAME INDIVIDUAL TRACERS
CC
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! State variables Info
        ! SubModel   Components  Name            Short Description:
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! pelagic  (O)              O2:   Oxygen
        ! pelagic  (P)              N1:   Phosphate
        ! pelagic  (N)              N3:   Nitrate
        ! pelagic  (N)              N4:   Ammonium
        ! pelagic  (Si)             N5:   Silicate

        integer ppO2o, ppN1p, ppN3n, ppN4n, ppN5s 
        parameter (ppO2o=1, ppN1p=2, ppN3n=3, ppN4n=4, ppN5s=5) 

