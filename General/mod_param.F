      MODULE modul_param

      IMPLICIT NONE

      public
CCC---------------------------------------------------------------------
CCC  OPA8, LODYC (15/11/96)
CCC---------------------------------------------------------------------
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
           
CC
CC	jpiglo  : first  dimension of global domain --> i	 
CC	jpjglo  : second dimension of global domain --> j
CC	jpk	: number of vertical levels
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
CC	jpi     : first  dimension of grid --> i	 
CC	jpi     : first  dimension of grid --> i	 
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
CC	jperio	: lateral cond. type for the global domain (4, 3, 2, 1
CC                or 0)  
CC	jpisl	: number of islands
CC	jpnisl	: maximum number of points per island 
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
CC	jpcpu	: maximum number of tasks
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
#    if defined key_trc_bfm
CCC---------------------------------------------------------------------
CCC
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
#    if defined key_INCLUDE_BFM_PELCO2

      INTEGER, parameter :: jptra = 51
      INTEGER, parameter :: jptra_dia = 10

#    else

      INTEGER, parameter :: jptra = 49
      INTEGER, parameter :: jptra_dia = 3
#    endif
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
        ! pelagic  (R)              N6:   Reduction Equivalents
        ! pelagic  (CNP)            B1:   Pelagic Bacteria
        ! pelagic  (CNPSiI)         P1:   Diatoms (group PhytoPlankton))
        ! pelagic  (CNPSiI)         P2:   Flagellates (group PhytoPlankton))
        ! pelagic  (CNPSiI)         P3:   PicoPhytoPlankton (group PhytoPlankton))
        ! pelagic  (CNPSiI)         P4:   Dinoflagellates (group PhytoPlankton))
        ! pelagic  (CNP)            Z3:   Carnivorous mesozooplankton (group MesoZooPlankton))
        ! pelagic  (CNP)            Z4:   Omnivorous mesozooplankton (group MesoZooPlankton))
        ! pelagic  (CNP)            Z5:   Microzooplankton (group MicroZooPlankton))
        ! pelagic  (CNP)            Z6:   Heterotrophic nanoflagellates (HNAN) (group MicroZooPlankton))
        ! pelagic  (CNPSi)          R1:   Labile Organic Carbon (LOC)
        ! pelagic  (C)              R2:   CarboHydrates (sugars)
        ! pelagic  (CNPSi)          R6:   Particulate Organic Carbon (POC)
        ! pelagic  (C)              R7:

C       INTEGER,parameter :: ppO2o=1, ppN1p=2, ppN3n=3, ppN4n=4, ppN5s=5, 
C    &     ppN6r=6, ppB1c=7, ppB1n=8, ppB1p=9, ppP1c=10, ppP2c=11, ppP3c=12, 
C    &     ppP4c=13, ppP1n=14, ppP2n=15, ppP3n=16, ppP4n=17, ppP1p=18, ppP2p=19, 
C    &     ppP3p=20, ppP4p=21, ppP1s=22, ppP2s=23, ppP3s=24, ppP4s=25, ppP1i=26, 
C    &     ppP2i=27, ppP3i=28, ppP4i=29, ppZ3c=30, ppZ4c=31, ppZ3n=32, ppZ4n=33, 
C    &     ppZ3p=34, ppZ4p=35, ppZ5c=36, ppZ6c=37, ppZ5n=38, ppZ6n=39, ppZ5p=40, 
C    &     ppZ6p=41, ppR1c=42, ppR1n=43, ppR1p=44, ppR1s=45, ppR2c=46, ppR6c=47, 
C    &     ppR6n=48, ppR6p=49, ppR6s=50, ppR7c=51

        integer,parameter :: ppO2o=1, ppN1p=2, ppN3n=3, ppN4n=4, ppO4n=5,
     &     ppN5s=6, ppN6r=7, ppB1c=8, ppB1n=9, ppB1p=10, ppP1c=11, ppP1n=12,
     &     ppP1p=13, ppP1i=14, ppP1s=15, ppP2c=16, ppP2n=17, ppP2p=18, ppP2i=19,
     &     ppP3c=20, ppP3n=21, ppP3p=22, ppP3i=23, ppP4c=24, ppP4n=25, ppP4p=26,
     &     ppP4i=27, ppZ3c=28, ppZ3n=29, ppZ3p=30, ppZ4c=31, ppZ4n=32, ppZ4p=33,
     &     ppZ5c=34, ppZ5n=35, ppZ5p=36, ppZ6c=37, ppZ6n=38, ppZ6p=39, ppR1c=40,
     &     ppR1n=41, ppR1p=42, ppR1s=43, ppR2c=44, ppR6c=45, ppR6n=46, ppR6p=47,
     &     ppR6s=48, ppR7c=49
C       diagnostic indexes
        integer,parameter :: ppPPG=1, ppNPP=2, ppBPP=3
#    if defined key_INCLUDE_BFM_PELCO2

        integer,parameter :: ppO3c=50, ppO3h=51
C       diagnostic indexes
        integer,parameter :: ppDIC=4, ppCO2=5,  ppHCO3=6, ppCO3=7, ppAc=8,
     &                       pppH=9, pppCO2=10

#    endif

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
