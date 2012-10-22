CC $Header: /cvsroot/opatm-bfm/opa_model/OPA/parameter.offline.h,v 1.2 2009-09-11 09:20:56 cvsogs01 Exp $
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
CC	original : 00 (O. Aumont M.-A. Foujols)
CC
CC Number of fluxes
CC ----------------
CC      jpflx   : number of fluxes
CC
      INTEGER jpflx
      PARAMETER(jpflx=13)
CC
CC	jptaux  : zonal wind stress
CC	jptauy  : meridional wind stress
CC	jpemp   : E - P in mm/day
CC	jpqsr   : solar radiation
CC
      INTEGER jpemp,jpqsr,jptaux,jptauy
      INTEGER jpwind,jpice
      PARAMETER(jptaux=1,jptauy=2,jpemp=4,jpqsr=6)
      PARAMETER(jpwind=3,jpice=5)
CC
CC
CC
CC Offline dynamic data size
CC -------------------------
CC      jpilec  : first horizontal dimension > or = jpi
CC      jpjlec  : second                     > or = jpj
CC      jpklec  : number of levels           > or = jpk
CC
      INTEGER jpilec,jpjlec,jpklec
      PARAMETER(jpilec=363,jpjlec=113,jpklec=32)
CC
