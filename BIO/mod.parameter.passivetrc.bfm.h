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
      INTEGER, parameter :: jptra = 51
CC      PARAMETER (jptra = 51)
CC
CC productive layer depth
CC ----------------------
CC       jpkb   : first vertical layers where biology is active
CC       jpkbm1 : jpkb - 1
CC
      INTEGER, parameter :: jpkb =72,
     &                      jpkbm1 =71
CC      PARAMETER (jpkb = 72,jpkbm1 = 71)
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

        integer,parameter :: ppO2o=1, ppN1p=2, ppN3n=3, ppN4n=4, ppN5s=5, 
     &     ppN6r=6, ppB1c=7, ppB1n=8, ppB1p=9, ppP1c=10, ppP2c=11, ppP3c=12, 
     &     ppP4c=13, ppP1n=14, ppP2n=15, ppP3n=16, ppP4n=17, ppP1p=18, ppP2p=19, 
     &     ppP3p=20, ppP4p=21, ppP1s=22, ppP2s=23, ppP3s=24, ppP4s=25, ppP1i=26, 
     &     ppP2i=27, ppP3i=28, ppP4i=29, ppZ3c=30, ppZ4c=31, ppZ3n=32, ppZ4n=33, 
     &     ppZ3p=34, ppZ4p=35, ppZ5c=36, ppZ6c=37, ppZ5n=38, ppZ6n=39, ppZ5p=40, 
     &     ppZ6p=41, ppR1c=42, ppR1n=43, ppR1p=44, ppR1s=45, ppR2c=46, ppR6c=47, 
     &     ppR6n=48, ppR6p=49, ppR6s=50, ppR7c=51

