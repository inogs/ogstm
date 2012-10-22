CCC---------------------------------------------------------------------
CCC
CCC                         PARAMETER.OFFLINE
CCC                       ********************
CCC
CCC  offline reading parameter
CCC
CCC  PURPOSE :
CCC  ---------
CCC     Include parameter FILE for offline specificity
CCC
CC   MODIFICATIONS :
CC   -------------
CC    original : 00 (O. Aumont M.-A. Foujols)
CC
CC Number of fluxes
CC ----------------
CC      jpflx   : number of fluxes
CC
      INTEGER :: jpflx 
CC
CC    jptaux  : zonal wind stress
CC    jptauy  : meridional wind stress
CC    jpemp   : E - P in mm/day
CC    jpqsr   : solar radiation
CC
      INTEGER :: jptaux,
     &           jptauy, 
     &           jpwind,
     &           jpemp,
     &           jpice,
     &           jpqsr
      
CC      PARAMETER(jptaux=1,jptauy=2,jpemp=4,jpqsr=6)
CC      PARAMETER(jpwind=3,jpice=5)
CC
CC
CC
CC Offline dynamic data size
CC -------------------------
CC      jpilec  : first horizontal dimension > or = jpi
CC      jpjlec  : second                     > or = jpj
CC      jpklec  : number of levels           > or = jpk
CC
      INTEGER :: jpilec,
     &           jpjlec,
     &           jpklec
CC      PARAMETER(jpilec=363,jpjlec=113,jpklec=32)
CC
