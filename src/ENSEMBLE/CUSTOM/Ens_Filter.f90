! developed by Simone Spada (sspada@ogs.it) at OGS

module Ens_Filter
    !use Seik, &
    use Seik_Likelihood, &
        only: Seik_Init, Seik_Finalize, Seik_Analysis
    
    implicit none
    
contains 
    
    subroutine Ens_Init_Filter
        call Seik_Init
    end subroutine
    
    subroutine Ens_Finalize_Filter
        call Seik_Finalize
    end subroutine
    
    subroutine EnsForecast
        
    end subroutine
    
    subroutine EnsAnalysis
        call Seik_Analysis
    end subroutine

end module
