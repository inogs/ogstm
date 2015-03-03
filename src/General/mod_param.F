      MODULE modul_param

      IMPLICIT NONE

      public
CC
CC Several domain sizes are parameterized
CC                                  global      domain   (jpiglo,jpjglo)
CC                                  computation domain   ( jpi  , jpj  )
CC
CC * if we dont use massively parallel computer (parameters jpni=jpnj=1)
CC      so jpiglo=jpi and jpjglo=jpj
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
CC
CC-CC      INTEGER,parameter :: jpni=2, 
CC-CC     &                     jpnj=1,   
CC-CC     &                     jpnij=jpni*jpnj,  
CC-CC     &                     jpreci=1, 
CC-CC     &                     jprecj=1  
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

      INTEGER :: jpiglo, jpjglo, jpk

CC      jpkmod  : modulo jpk and jpnij
CC
CC-CC      INTEGER,parameter :: jpkmod=(1+(jpk-1)/jpnij)
      INTEGER :: jpkmod

CC Matrix size
CC -----------
CC    jpi     : first  dimension of grid --> i
CC    jpi     : first  dimension of grid --> i
CC
CC-CC      INTEGER,parameter :: jpi=(jpiglo-2*jpreci + (jpni-1))/jpni + 2*jpreci
CC-CC      INTEGER,parameter :: jpj=(jpjglo-2*jprecj + (jpnj-1))/jpnj + 2*jprecj
      INTEGER :: jpi
      INTEGER :: jpj
CC
CC
CC Other dimension parameters
CC --------------------------
CC      jpim1   :  jpi - 1
CC      jpjm1   :  jpj - 1
CC      jpkm1   :  jpk - 1
CC      jpij    :  jpi x jpj
CC
CC-CC      INTEGER,parameter :: jpim1=jpi-1,  
CC-CC     &                     jpjm1=jpj-1,  
CC-CC     &                     jpkm1=jpk-1,  
CC-CC     &                     jpij=jpi*jpj 
      INTEGER :: jpim1,  
     &           jpjm1,  
     &           jpkm1,  
     &           jpij 
CC

CC
CC
CC Domain characteristics
CC ----------------------
CC    jperio    : lateral cond. type for the global domain (4, 3, 2, 1
CC                or 0)  
CC
CC
CC
         INTEGER :: jperio  


CC
CC Global grid characteristics
CC ---------------------------
CC
CC Offline reading parameter
CC -------------------------
CC
!#include "mod.param.off.h"
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



!end include "mod.param.off.h"
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
CC      original : 00 (O. Aumont M.-A. Foujols)
CC
CC Number of fluxes
CC ----------------
CC      jpflx   : number of fluxes
CC
CC
CC      jptaux  : zonal wind stress
CC      jptauy  : meridional wind stress
CC      jpemp   : E - P in mm/day
CC      jpqsr   : solar radiation
       INTEGER :: jpkef

CC Passive tracers parameter
CC -------------------------
CC
#ifdef key_passivetrc
#    if defined key_trc_bfm


#include "BFM_var_list.h"
CC
CC productive layer depth
CC ----------------------
CC       jpkb   : first vertical layers where biology is active
CC       jpkbm1 : jpkb - 1
CC
      INTEGER :: jpkb ,
     &           jpkbm1
CC
CC number of biological trends
CC ---------------------------
CC
      INTEGER, parameter :: jpdiabio = -1


#    elif defined key_test00
      INTEGER :: jpkb ,jpkbm1
#    else
CC    default CASE : temperature and salinity as passive tracers
      INTEGER, parameter :: jptra = 2
CC      PARAMETER (jptra = 2)
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
      INTEGER, parameter :: jptra = 5
CC
CC no passive tracer 
CC
#endif
CC

         END MODULE 
