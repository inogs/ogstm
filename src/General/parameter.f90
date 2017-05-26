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

#include "BFM_var_list.h"

!! productive layer depth

      INTEGER :: jpkb , & ! : first vertical layers where biology is active
     &           jpkbm1 !  jpkb - 1

#endif


      END MODULE modul_param
