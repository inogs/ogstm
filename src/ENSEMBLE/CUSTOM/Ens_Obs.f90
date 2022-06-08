! developed by Simone Spada (sspada@ogs.it) at OGS

module Ens_Obs
    use Ens_Custom, &   
        only: nk_DAstate, nj_DAstate, ni_DAstate, ntra_DAstate
    use ObsSatellite, &
        only: Sat_Init, Sat_Finalize, Sat_LoadObs, Sat_Misfit, Sat_H, Sat_SqrtR1
    
    implicit none
    
    logical, parameter :: UseSatellite=.true.
    
    integer :: nObs, n_SqrtR1
    integer :: Sat_nObs, Sat_n_SqrtR1

contains

    Subroutine Ens_Init_Obs
        
        call Sat_Init
        
    end Subroutine
    
    subroutine Ens_Finalize_Obs
        use mpi
        
        integer ierror
        
        call sat_Finalize
        
    end subroutine
    
    subroutine Ens_LoadObs(DateString)
        
        character(len=*), intent(in) :: DateString
        
        nObs=0 
        n_SqrtR1=0
        Sat_nObs=0
        Sat_n_SqrtR1=0
        
        if (UseSatellite) then
            call Sat_LoadObs(DateString, Sat_nObs, Sat_n_SqrtR1)
            nObs=nObs+Sat_nObs
            n_SqrtR1=n_SqrtR1+Sat_n_SqrtR1
        end if
        
    end subroutine
    
    function Misfit(ObsState)
    
        double precision, dimension(nObs,nj_DAstate, ni_DAstate), intent(in) :: ObsState
        double precision, dimension(nObs,nj_DAstate, ni_DAstate) :: Misfit
        
        Misfit(1:Sat_nObs,:,:)=Sat_Misfit(ObsState(1:Sat_nObs,:,:))
        
    end function
    
    function H(State)
        
        double precision, dimension(nk_DAstate,nj_DAstate,ni_DAstate,ntra_DAstate), intent(in) :: State
        double precision, dimension(nObs,nj_DAstate,ni_DAstate) :: H
        
        H(1:Sat_nObs,:,:)=Sat_H(State)
        
    end function
    
    function SqrtR1(HLi)
        
        double precision, dimension(nObs, nj_DAstate, ni_DAstate), intent(in) :: HLi
        double precision, dimension(n_SqrtR1, nj_DAstate, ni_DAstate) :: SqrtR1
        
        SqrtR1(1:Sat_n_SqrtR1,:,:)=Sat_SqrtR1(HLi(1:Sat_nObs,:,:))
        
    end function

end module

