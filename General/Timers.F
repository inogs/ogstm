      MODULE Timers

CC----------------------------------------------------------------------
CC      cronometer for function calls
CC -----------------------------------

!      REAL(8) :: dayparttime         =0.0,  daytottime       = 0.0
      REAL(8) :: forcing_phys_partTime=0.0, forcing_phys_TotTime = 0.0
      REAL(8) :: forcing_kext_partTime=0.0, forcing_kext_TotTime = 0.0
      REAL(8) :: bc_Co2_partTime      =0.0, bc_Co2_TotTime       = 0.0
      REAL(8) :: bc_atm_partTime      =0.0, bc_atm_TotTime       = 0.0
      REAL(8) :: bc_tin_partTime      =0.0, bc_tin_TotTime       = 0.0
      REAL(8) :: bc_gib_partTime      =0.0, bc_gib_TotTime       = 0.0
      REAL(8) :: density_partTime     =0.0, density_TotTime      = 0.0
      REAL(8) :: ave_partTime         =0.0, ave_TotTime          = 0.0
      REAL(8) :: flx_partTime         =0.0, flx_TotTime          = 0.0
      REAL(8) :: trcdiaparttime       =0.0, trcdiatottime     = 0.0
      REAL(8) :: trcstpparttime       =0.0, trcstptottime     = 0.0
      REAL(8) :: trcwriparttime       =0.0, trcwritottime     = 0.0
      REAL(8) :: trcsmsparttime       =0.0, trcsmstottime     = 0.0
      REAL(8) :: trcadvparttime       =0.0, trcadvtottime     = 0.0
      REAL(8) :: trcdmpparttime       =0.0, trcdmptottime     = 0.0
      REAL(8) :: trclaphdfparttime    =0.0, trclaphdftottime  = 0.0
      REAL(8) :: trcbilaphdfparttime  =0.0, trcbilaphdftottime= 0.0
      REAL(8) :: trcbblparttime       =0.0, trcbbltottime     = 0.0
      REAL(8) :: trczdfparttime       =0.0, trczdftottime     = 0.0
      REAL(8) :: trcnxtparttime       =0.0, trcnxttottime     = 0.0
      REAL(8) :: trcoptparttime       =0.0, trcopttottime     = 0.0
      REAL(8) :: trcbioparttime       =0.0, trcbiotottime     = 0.0
      REAL(8) :: stpparttime          =0.0, stptottime        = 0.0
      REAL(8) :: BIOparttime          =0.0, BIOtottime        = 0.0
      REAL(8) :: snutelparttime       =0.0, snuteltottime     = 0.0
      REAL(8) :: trcsbcparttime       =0.0, trcsbctottime     = 0.0
      CONTAINS





      SUBROUTINE reset_Timers()
      forcing_phys_TotTime = 0.0
      forcing_kext_TotTime = 0.0
      bc_Co2_TotTime       = 0.0
      bc_atm_TotTime       = 0.0
      bc_tin_TotTime       = 0.0
      bc_gib_TotTime       = 0.0
      density_TotTime      = 0.0
      ave_TotTime          = 0.0
      flx_TotTime          = 0.0
      trcdiatottime        = 0.0
      trcstptottime        = 0.0
      trcwritottime        = 0.0
      trcsmstottime        = 0.0
      trcadvtottime        = 0.0
      trcdmptottime        = 0.0
      trclaphdftottime     = 0.0
      trcbilaphdftottime   = 0.0
      trcbbltottime        = 0.0
      trczdftottime        = 0.0
      trcnxttottime        = 0.0
      trcopttottime        = 0.0
      trcbiotottime        = 0.0
      stptottime           = 0.0
      BIOtottime           = 0.0
      snuteltottime        = 0.0
      trcsbctottime        = 0.0
      END SUBROUTINE reset_Timers

      END MODULE TIMERS



