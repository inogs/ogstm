! developed by Simone Spada (sspada@ogs.it) at OGS

module Ens_Params
    use global_mem, &
        only: RLEN
    use mem, &
        only: iiB1, iiZ5, iiZ6, iiP1, iiP3, iiR6
    use mem_PelBac, &
        only: p_pu_ra
    use mem_MicroZoo, &
        only: p_sum
    use mem_Phyto, &
        only: p_qlcPPY, p_qpcPPY, p_qplc, p_srs
    use myalloc, &
        only: vsed
        
    implicit none
    
    integer, parameter :: n_params=9
    
    real(RLEN), dimension(:), allocatable :: ParamArray
    character(len=100), dimension(n_params) :: ParamNames
    
contains
    
    subroutine Ens_Init_Params
        allocate(ParamArray(n_params))
        
        ParamNames(1)='PRMS.B1.p_pu_ra'
        ParamNames(2)='PRMS.Z5.p_sum'
        ParamNames(3)='PRMS.Z6.p_sum'
        ParamNames(4)='PRMS.P1.p_qlcPPY'
        ParamNames(5)='PRMS.P3.p_qlcPPY'
        ParamNames(6)='PRMS.P3.p_qpcPPY'
        ParamNames(7)='PRMS.P1.p_qplc'
        ParamNames(8)='PRMS.P1.p_srs'
        ParamNames(9)='PRMS.R6.rm'
        
    end subroutine

    subroutine Ens_Finalize_Params
        deallocate(ParamArray)
        
    end subroutine
    
    subroutine Ens_SetParams
        p_pu_ra(iiB1)=ParamArray(1)
        p_sum(iiZ5)=ParamArray(2)
        p_sum(iiZ6)=ParamArray(3)
        p_qlcPPY(iiP1)=ParamArray(4)
        p_qlcPPY(iiP3)=ParamArray(5)
        p_qpcPPY(iiP3)=ParamArray(6)
        p_qplc(iiP1)=ParamArray(7)
        p_srs(iiP1)=ParamArray(8)
        vsed=ParamArray(9)
        
    end subroutine


end module
