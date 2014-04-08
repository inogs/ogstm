      MODULE modul_param

      IMPLICIT NONE

      public
CC
CC Several domain sizes are parameterized
CC                                  global      domain   (jpiglo,jpjglo)
CC                                  computation domain   ( jpi  , jpj  )
CC                                  data        domain   (jpidta,jpjdta)
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
CC      jpidta  : first horizontal dimension > or = jpi
CC      jpjdta  : second                     > or = jpj
CC      jpkdta  : number of levels           > or = jpk
CC

CCC 25/2/2004 Paolo  Dimension of Mediterranean Grid
      INTEGER :: jpiglo,    
     &           jpjglo,    
     &           jpk 

      INTEGER :: jpidta,  
     &           jpjdta, 
     &           jpkdta
CCC
CC
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
CC    jpisl    : number of islands
CC    jpnisl    : maximum number of points per island
CC
CC
CC
         INTEGER :: jperio,  
     &              jpisl,   
     &              jpnisl 

CC
CC
CC Multitasking parameter
CC ----------------------
CC    jpcpu    : maximum number of tasks
CC      jpbyt   : number of bytes for floating  4 or 8 
CC      jpbytda : number of bytes for floating in data files 4 or 8
CC
      INTEGER :: jpcpu=8,  
     &           jpbyt=8,  
     &           jpbi3e=4, 
     &           jpbytda=4

CC
CC Global grid characteristics
CC ---------------------------
CC
CC Offline reading parameter
CC -------------------------
CC
#include "mod.param.off.h"
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
CC
CC
CC
       INTEGER :: jpkef
CC
CC Offline dynamic data size
CC -------------------------
CC      jpilec  : first horizontal dimension > or = jpi
CC      jpjlec  : second                     > or = jpj
CC      jpklec  : number of levels           > or = jpk
CC
CC
CC Passive tracers parameter
CC -------------------------
CC
#ifdef key_passivetrc
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
CC    original : 96 (M. Levy)
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
#    if defined key_trc_bfm
CCC---------------------------------------------------------------------
CCC                         PARAMETER passivetrc.bfm
CCC                       *****************************
CCC
CCC  purpose :
CCC  ---------
CCC     INCLUDE PARAMETER FILE for passive tracer BFM model
CCC
CCC  --------------
CC
CC       jptra  : number of tracers
CC

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
CC      PARAMETER (jpdiabio = -1)
CC

#    elif defined key_test00
CC    test CASE : nutrient adv-dif without reaction
      INTEGER, parameter :: jptra = 5
      INTEGER, parameter :: ppO2o=1, ppN1p=2, ppN3n=3, ppN4n=4, ppN5s=5
CC productive layer depth
CC ----------------------
CC       jpkb   : first vertical layers where biology is active
CC       jpkbm1 : jpkb - 1
CC
      INTEGER :: jpkb ,
     &           jpkbm1

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
