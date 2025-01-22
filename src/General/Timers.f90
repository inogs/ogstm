       MODULE Timers

!!----------------------------------------------------------------------
!!      cronometer for function calls
!! -----------------------------------


      double precision :: forcing_phys_partTime=0.0, forcing_phys_TotTime = 0.0
      double precision :: forcing_kext_partTime=0.0, forcing_kext_TotTime = 0.0
      double precision :: bc_Co2_partTime      =0.0, bc_Co2_TotTime       = 0.0
      double precision :: bc_atm_Hg0_partTime  =0.0, bc_atm_Hg0_TotTime   = 0.0
      double precision :: bc_atm_partTime      =0.0, bc_atm_TotTime       = 0.0
      double precision :: bc_tin_partTime      =0.0, bc_tin_TotTime       = 0.0
      double precision :: bc_gib_partTime      =0.0, bc_gib_TotTime       = 0.0
      double precision :: bc_dar_partTime      =0.0, bc_dar_TotTime       = 0.0
      double precision :: density_partTime     =0.0, density_TotTime      = 0.0
      double precision :: ave_partTime         =0.0, ave_TotTime          = 0.0
      double precision :: flx_partTime         =0.0, flx_TotTime          = 0.0
      double precision :: trcdiaparttime       =0.0, trcdiatottime     = 0.0
      double precision :: trcstpparttime       =0.0, trcstptottime     = 0.0
      double precision :: trcwriparttime       =0.0, trcwritottime     = 0.0
      double precision :: trcsmsparttime       =0.0, trcsmstottime     = 0.0
      double precision :: trcadvparttime       =0.0, trcadvtottime     = 0.0
      double precision :: trcdmpparttime       =0.0, trcdmptottime     = 0.0
      double precision :: trclaphdfparttime    =0.0, trclaphdftottime  = 0.0
      double precision :: trcbilaphdfparttime  =0.0, trcbilaphdftottime= 0.0
      double precision :: trczdfparttime       =0.0, trczdftottime     = 0.0
      double precision :: trcnxtparttime       =0.0, trcnxttottime     = 0.0
      double precision :: trcoptparttime       =0.0, trcopttottime     = 0.0
      double precision :: trcbioparttime       =0.0, trcbiotottime     = 0.0
      double precision :: stpparttime          =0.0, stptottime        = 0.0
      double precision :: BIOparttime          =0.0, BIOtottime        = 0.0
      double precision :: snutelparttime       =0.0, snuteltottime     = 0.0
      double precision :: trcsbcparttime       =0.0, trcsbctottime     = 0.0
      double precision :: checkVparttime       =0.0, checkVtottime     = 0.0
      double precision :: DAparttime           =0.0, DAtottime         = 0.0
      CONTAINS





      SUBROUTINE reset_Timers()
      forcing_phys_TotTime = 0.0
      forcing_kext_TotTime = 0.0
      bc_Co2_TotTime       = 0.0
      bc_atm_Hg0_TotTime   = 0.0
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
      trczdftottime        = 0.0
      trcnxttottime        = 0.0
      trcopttottime        = 0.0
      trcbiotottime        = 0.0
      stptottime           = 0.0
      BIOtottime           = 0.0
      snuteltottime        = 0.0
      trcsbctottime        = 0.0
      checkVtottime        = 0.0
      END SUBROUTINE reset_Timers

      END MODULE TIMERS



