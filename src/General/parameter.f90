      MODULE modul_param

      !use mem

      IMPLICIT NONE

      public

!! Several domain sizes are parameterized
!!                                  global      domain   (jpjglo,jpiglo)
!!                                  computation domain   ( jpi  , jpj  )
!!
!!
!! Large domain matrix size
!! ------------------------
!!      jpni    : number of processors following i 
!!      jpnj    : number of processors following j
!!      jpnij   : < or = jpni x jpnj number of processors 
!!                (i.e. local domains)
!!                One can avoid land processors
!!      jpreci  : number of lines for overlap 
!!      jprecj  : number of lines for overlap 

      INTEGER :: jpni, &
     &           jpnj, &
     &           jpnij, &
     &           jpreci, &
     &           jprecj 
           
!!    jpiglo  : first  dimension of global domain --> i
!!    jpjglo  : second dimension of global domain --> j
!!    jpk    : number of vertical levels
!!
!! Original data size
!! ------------------

      INTEGER :: jpk, jpjglo, jpiglo


!! Matrix size
      INTEGER :: jpi,jpj
      INTEGER :: jpim1,  jpjm1,  jpkm1,  jpij
!! Domain characteristics
!! ----------------------

      INTEGER :: jperio
      INTEGER :: jpflx ! number of fluxes

!!    jpemp   : E - P in mm/day
!!    jpqsr   : solar radiation
      INTEGER :: jpwind,jpemp,jpice,jpqsr,jpkef


!! Passive tracers parameter
#ifdef key_trc_bfm

      !! WARNING the var below must be become input parameter
!      INTEGER, parameter :: jptra = 51
!      INTEGER, parameter :: jptra_dia = 21
!      INTEGER, parameter :: jptra_dia_2d = 1

!#include "BFM_var_list.h"


!! productive layer depth

      INTEGER :: jpkb , & ! : first vertical layers where biology is active
     &           jpkbm1 !  jpkb - 1

#endif




      INTEGER, parameter :: jptra = 1 !51

      INTEGER, parameter :: jptra_var = 100

      INTEGER, parameter :: jptra_flux = 2

      INTEGER, parameter :: jptra_dia = 0 !jptra_var + jptra_flux
       
      INTEGER, parameter :: jptra_dia_2d = 0 !11

 
! State variables

        integer,parameter :: ppO2o=1, ppN1p=2, ppN3n=3, ppN4n=4, ppO4n=5
     integer,parameter ::  ppN5s=6, ppN6r=7, ppB1c=8, ppB1n=9, ppB1p=10, ppP1c=11, ppP1n=12, ppP1p=13
     integer,parameter ::  ppP1l=14, ppP1s=15, ppP2c=16, ppP2n=17, ppP2p=18, ppP2l=19, ppP2s=0
     integer,parameter ::  ppP3c=20, ppP3n=21, ppP3p=22, ppP3l=23, ppP3s=0, ppP4c=24, ppP4n=25
     integer,parameter ::  ppP4p=26, ppP4l=27, ppP4s=0, ppZ3c=28, ppZ3n=29, ppZ3p=30, ppZ4c=31
     integer,parameter ::  ppZ4n=32, ppZ4p=33, ppZ5c=34, ppZ5n=35, ppZ5p=36, ppZ6c=37, ppZ6n=38
     integer,parameter ::  ppZ6p=39, ppR1c=40, ppR1n=41, ppR1p=42, ppR1s=0, ppR2c=43, ppR2n=0
     integer,parameter ::  ppR2p=0, ppR2s=0, ppR3c=44, ppR3n=0, ppR3p=0, ppR3s=0, ppR6c=45, ppR6n=46
     integer,parameter ::  ppR6p=47, ppR6s=48, ppO3c=49, ppO3h=50, ppO5c=51, ppO5h=0


!       diagnostic indexes
        integer,parameter :: ppETW=1, ppESW=2, ppERHO=3, ppEIR=4, ppESS=5
     integer,parameter ::  ppEPR=6, ppDepth=7, ppVolume=8, ppArea=9, ppDIC=10, ppCO2=11, pppCO2=12
     integer,parameter ::  ppHCO3=13, ppCO3=14, ppALK=15, pppH=16, ppOCalc=17, ppOArag=18
     integer,parameter ::  pptotpelc=19, pptotpeln=20, pptotpelp=21, pptotpels=22, ppcxoO2=23
     integer,parameter ::  ppeO2mO2=24, ppChla=25, ppffCO2=26, ppflPTN6r=27, ppflN3O4n=28
     integer,parameter ::  ppflN4N3n=29, ppsediR2=30, ppsediR6=31, ppsediO5=32, ppxEPS=33
     integer,parameter ::  ppABIO_eps=34, ppqpcPPY_iiP1=35, ppqpcPPY_iiP2=36, ppqpcPPY_iiP3=37
     integer,parameter ::  ppqpcPPY_iiP4=38, ppqncPPY_iiP1=39, ppqncPPY_iiP2=40, ppqncPPY_iiP3=41
     integer,parameter ::  ppqncPPY_iiP4=42, ppqscPPY_iiP1=43, ppqscPPY_iiP2=44, ppqscPPY_iiP3=45
     integer,parameter ::  ppqscPPY_iiP4=46, ppqlcPPY_iiP1=47, ppqlcPPY_iiP2=48, ppqlcPPY_iiP3=49
     integer,parameter ::  ppqlcPPY_iiP4=50, ppqccPPY_iiP1=51, ppqccPPY_iiP2=52, ppqccPPY_iiP3=53
     integer,parameter ::  ppqccPPY_iiP4=54, ppBFM1D_exR2ac_iiP1=55, ppBFM1D_exR2ac_iiP2=56
     integer,parameter ::  ppBFM1D_exR2ac_iiP3=57, ppBFM1D_exR2ac_iiP4=58, ppqpcMEZ_iiZ3=59
     integer,parameter ::  ppqpcMEZ_iiZ4=60, ppqncMEZ_iiZ3=61, ppqncMEZ_iiZ4=62, ppqpcMIZ_iiZ5=63
     integer,parameter ::  ppqpcMIZ_iiZ6=64, ppqncMIZ_iiZ5=65, ppqncMIZ_iiZ6=66, ppqpcOMT_iiR1=67
     integer,parameter ::  ppqpcOMT_iiR2=68, ppqpcOMT_iiR3=69, ppqpcOMT_iiR6=70, ppqncOMT_iiR1=71
     integer,parameter ::  ppqncOMT_iiR2=72, ppqncOMT_iiR3=73, ppqncOMT_iiR6=74, ppqscOMT_iiR1=75
     integer,parameter ::  ppqscOMT_iiR2=76, ppqscOMT_iiR3=77, ppqscOMT_iiR6=78, ppqpcPBA_iiB1=79
     integer,parameter ::  ppqncPBA_iiB1=80, ppsediPPY_iiP1=81, ppsediPPY_iiP2=82, ppsediPPY_iiP3=83
     integer,parameter ::  ppsediPPY_iiP4=84, ppsediMIZ_iiZ5=85, ppsediMIZ_iiZ6=86, ppsediMEZ_iiZ3=87
     integer,parameter ::  ppsediMEZ_iiZ4=88, ppsunPPY_iiP1=89, ppsunPPY_iiP2=90, ppsunPPY_iiP3=91
     integer,parameter ::  ppsunPPY_iiP4=92, ppeiPPY_iiP1=93, ppeiPPY_iiP2=94, ppeiPPY_iiP3=95
     integer,parameter ::  ppeiPPY_iiP4=96, ppELiPPY_iiP1=97, ppELiPPY_iiP2=98, ppELiPPY_iiP3=99
     integer,parameter ::  ppELiPPY_iiP4=100


!       flux indexes
        integer, parameter:: ppruPPYc=1, ppresPPYc=2


!       variables 2d
        integer, parameter:: ppEPCO2air=1, ppCO2airflux=2, ppArea2d=3
     integer,parameter ::  ppThereIsLight=4, ppSUNQ=5, ppEWIND=6, pptotsysc=7, pptotsysn=8
     integer,parameter ::  pptotsysp=9, pptotsyss=10, ppEICE=11


      END MODULE modul_param
