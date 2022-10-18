! developed by Simone Spada (sspada@ogs.it) at OGS

module Ens_Params
    !use global_mem, &
    !    only: RLEN
    use mem, &
        only: iiB1, iiZ5, iiZ6, iiP1, iiP3, iiR6
    use mem_PelBac, &
        only: p_pu_ra
    use mem_MicroZoo, &
        only: p_sum
    use mem_Phyto, &
        only: p_qlcPPY, p_qpcPPY, p_qplc, p_srs
    use myalloc, &
        only: lwp, &
            vsed
        
    use Ens_Mem, &
        only: EnsIOUnit
        
    implicit none
    
    integer, parameter :: n_params=9
    
    double precision, dimension(:), allocatable :: ParamArray
    
    character(len=100), dimension(:), allocatable :: ParamName
    logical, dimension(:), allocatable :: UseParam
    
contains
    
    subroutine Ens_Init_Params
    
        integer indexi
    
        namelist/Params_setup/ ParamName, UseParam
        
        allocate(ParamArray(n_params))
        allocate(ParamName(n_params))
        allocate(UseParam(n_params))
        
        OPEN(unit=EnsIOUnit, file='namelist.params', status='OLD')
        
            ParamName(1)='PRMS.B1.p_pu_ra'
            UseParam(1)=.true.
            
            ParamName(2)='PRMS.Z5.p_sum'
            UseParam(2)=.true.
            
            ParamName(3)='PRMS.Z6.p_sum'
            UseParam(3)=.true.
            
            ParamName(4)='PRMS.P1.p_qlcPPY'
            UseParam(4)=.false.
            
            ParamName(5)='PRMS.P3.p_qlcPPY'
            UseParam(5)=.false.
            
            ParamName(6)='PRMS.P3.p_qpcPPY'
            UseParam(6)=.true.
            
            ParamName(7)='PRMS.P1.p_qplc'
            UseParam(7)=.true.
            
            ParamName(8)='PRMS.P1.p_srs'
            UseParam(8)=.true.
            
            ParamName(9)='PRMS.R6.rm'
            UseParam(9)=.true.
            
            REWIND(EnsIOUnit)
            READ(EnsIOUnit, Params_setup)
            
            IF(lwp) THEN
                WRITE(*,*) ''
                WRITE(*,*) 'Params_setup'
                WRITE(*,*) ''
                
                do indexi=1, n_params                
                    WRITE(*,*) ' ParamName(',indexi,'): ', trim(ParamName(indexi))
                    WRITE(*,*) ' UseParam(',indexi,'): ', UseParam(indexi)
                    WRITE(*,*) ''                
                end do
                
                WRITE(*,*) ''
            END IF

        CLOSE(EnsIOUnit)
        
    end subroutine

    subroutine Ens_Finalize_Params
    
        deallocate(ParamArray)
        deallocate(ParamName)
        deallocate(UseParam)
        
    end subroutine
    
    subroutine Ens_SetParams
    
        if (UseParam(1)) p_pu_ra(iiB1)=ParamArray(1)
        if (UseParam(2)) p_sum(iiZ5)=ParamArray(2)
        if (UseParam(3)) p_sum(iiZ6)=ParamArray(3)
        if (UseParam(4)) p_qlcPPY(iiP1)=ParamArray(4)
        if (UseParam(5)) p_qlcPPY(iiP3)=ParamArray(5)
        if (UseParam(6)) p_qpcPPY(iiP3)=ParamArray(6)
        if (UseParam(7)) p_qplc(iiP1)=ParamArray(7)
        if (UseParam(8)) p_srs(iiP1)=ParamArray(8)
        if (UseParam(9)) vsed=ParamArray(9)
        
    end subroutine
      
    subroutine Ens_ReadParams(prefix)
        !use Ens_Mem, &
        !    only: EnsIOUnit
        use TIME_MANAGER, &
            only: DateStart
        !use Ens_Params, &
        !    only: n_params, &
        !        ParamArray, ParamName
                
        character(LEN=*), intent(in) :: prefix
        
        integer indexi
        
        do indexi=1,n_params
            if (UseParam(indexi)) then
                open(EnsIOUnit, file=trim(prefix)//'.'//DateStart//'.'//trim(ParamName(indexi))//'.txt', status = 'old')
                    read(EnsIOUnit,*) ParamArray(indexi)
                close(EnsIOUnit)
            end if
        end do
        
    end subroutine


end module
