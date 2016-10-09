      MODULE modul_param

      IMPLICIT NONE

      public

CC Several domain sizes are parameterized
CC                                  global      domain   (jpjglo,jpiglo)
CC                                  computation domain   ( jpi  , jpj  )
CC
CC
CC Large domain matrix size
CC ------------------------
CC      jpni    : number of processors following i 
CC      jpnj    : number of processors following j
CC      jpnij   : < or = jpni x jpnj number of processors 
CC                (i.e. local domains)
CC                One can avoid land processors
CC      jpreci  : number of lines for overlap 
CC      jprecj  : number of lines for overlap 

      INTEGER :: jpni,
     &           jpnj,
     &           jpnij,
     &           jpreci,
     &           jprecj
           
CC    jpiglo  : first  dimension of global domain --> i
CC    jpjglo  : second dimension of global domain --> j
CC    jpk    : number of vertical levels
CC
CC Original data size
CC ------------------

      INTEGER :: jpk, jpjglo, jpiglo


CC Matrix size
      INTEGER :: jpi,jpj
      INTEGER :: jpim1,  jpjm1,  jpkm1,  jpij
CC Domain characteristics
CC ----------------------

      INTEGER :: jperio
      INTEGER :: jpflx ! number of fluxes

CC    jpemp   : E - P in mm/day
CC    jpqsr   : solar radiation
      INTEGER :: jpwind,jpemp,jpice,jpqsr,jpkef


CC Passive tracers parameter
#if defined key_trc_bfm


#include "BFM_var_list.h"

CC productive layer depth

      INTEGER :: jpkb , ! : first vertical layers where biology is active
     &           jpkbm1 !  jpkb - 1

#endif


      END MODULE modul_param
