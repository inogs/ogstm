! developed by Simone Spada (sspada@ogs.it) at OGS

module Ens_Params
    use mpi
    
    !use global_mem, &
    !    only: RLEN
    use modul_param, &
        only: jpi, jpj
    use mem, &
        only: iiB1, iiZ5, iiZ6, iiP1, iiP3, iiR6
    use mem_PelBac, &
        only: p_pu_ra
    use mem_MicroZoo, &
        only: p_sum
    use mem_Phyto, &
        only: p_qlcPPY, p_qpcPPY, p_qplc, p_srs
    use myalloc, &
        only: lwp, nldi, nldj, nlei, nlej, &
            vsed
            
    use ObsSatellite_Likelihood, &
        only: SatMultError, SatAddError, SatBias, &
            Sat_std_additive, Sat_std_log, Sat_bias, win_Input_Sat
    use Ens_Custom, &
        only: Inflation, win_Input_Custom
        
    use Ens_Mem, &
        only: EnsSize, EnsRankZero, &
            ForgettingFactor, &
            EnsIOUnit, &
            Ens_shared_alloc
    use ogstm_mpi_module, &
        only: glcomm
    !use IO_module, only: readnc_slice_double_2d !!! module do not exists, but it should be created.
        
    implicit none
    
    integer, parameter :: n_params=13
    
!     double precision, dimension(:,:,:), allocatable, target :: ParamArray
    
    double precision, pointer, contiguous, dimension(:,:,:) :: ParamArray
    integer :: n_ParamArray
    integer :: win_ParamArray
    double precision, pointer, contiguous, dimension(:,:,:,:) :: gl_ParamArray
    
    
    character(len=100), dimension(:), allocatable :: ParamName
    logical, dimension(:), allocatable :: UseParam, SaveParam
    integer, dimension(:), allocatable :: ParamSubset
    
contains
    
subroutine Ens_Init_Params
    
    double precision, dimension(:), POINTER, contiguous :: member_pointer
    double precision, dimension(:,:), POINTER, contiguous :: global_pointer
    
    integer indexi

    namelist/Params_setup/ ParamName, UseParam, ParamSubset, SaveParam
    
    n_ParamArray=jpj*jpi*n_params
    call Ens_shared_alloc(n_ParamArray, member_pointer, global_pointer, win_ParamArray)
    ParamArray(1:n_params, 1:jpj, 1:jpi)=>member_pointer
    gl_ParamArray(1:n_params, 1:jpj, 1:jpi, 0:EnsSize-1)=>global_pointer
    ParamArray  = huge(ParamArray(1,1,1))
    
!     allocate(ParamArray(n_params, jpj, jpi))
    allocate(ParamName(n_params))
    allocate(UseParam(n_params))
    allocate(ParamSubset(n_params))
    allocate(SaveParam(n_params))
    
    OPEN(unit=EnsIOUnit, file='namelist.params', status='OLD')
    
        ParamName(1)='PRMS.B1.p_pu_ra'
        UseParam(1)=.true.
        ParamSubset(1)=0
        SaveParam(1)=.false.
        
        ParamName(2)='PRMS.Z5.p_sum'
        UseParam(2)=.true.
        ParamSubset(2)=0
        SaveParam(2)=.false.
        
        ParamName(3)='PRMS.Z6.p_sum'
        UseParam(3)=.true.
        ParamSubset(3)=0
        SaveParam(3)=.false.
        
        ParamName(4)='PRMS.P1.p_qlcPPY'
        UseParam(4)=.true.
        ParamSubset(4)=0
        SaveParam(4)=.false.
        
        ParamName(5)='PRMS.P3.p_qlcPPY'
        UseParam(5)=.true.
        ParamSubset(5)=0
        SaveParam(5)=.false.
        
        ParamName(6)='PRMS.P3.p_qpcPPY'
        UseParam(6)=.true.
        ParamSubset(6)=0
        SaveParam(6)=.false.
        
        ParamName(7)='PRMS.P1.p_qplc'
        UseParam(7)=.true.
        ParamSubset(7)=0
        SaveParam(7)=.false.
        
        ParamName(8)='PRMS.P1.p_srs'
        UseParam(8)=.true.
        ParamSubset(8)=0
        SaveParam(8)=.false.
        
        ParamName(9)='PRMS.R6.rm'
        UseParam(9)=.true.
        ParamSubset(9)=0
        SaveParam(9)=.false.
        
        ParamName(10)='PRMS.Infl'
        UseParam(10)=.true.
        ParamSubset(10)=1
        SaveParam(10)=.true.
        
        ParamName(11)='PRMS.Sat.Mult'
        UseParam(11)=.true.
        ParamSubset(11)=1
        SaveParam(11)=.true.
        
        ParamName(12)='PRMS.Sat.Add'
        UseParam(12)=.true.
        ParamSubset(12)=1
        SaveParam(12)=.true.
        
        ParamName(13)='PRMS.Sat.Bias'
        UseParam(13)=.false.
        ParamSubset(13)=1
        SaveParam(13)=.true.
        
        REWIND(EnsIOUnit)
        READ(EnsIOUnit, Params_setup)
        
        IF(lwp) THEN
            WRITE(*,*) ''
            WRITE(*,*) 'Params_setup'
            WRITE(*,*) ''
            
            do indexi=1, n_params                
                WRITE(*,*) ' ParamName(',indexi,'): ', trim(ParamName(indexi))
                WRITE(*,*) ' UseParam(',indexi,'): ', UseParam(indexi)
                WRITE(*,*) ' ParamSubset(',indexi,'): ', ParamSubset(indexi)
                WRITE(*,*) ' SaveParam(',indexi,'): ', SaveParam(indexi)
                WRITE(*,*) ''                
            end do
            
            WRITE(*,*) ''
        END IF

    CLOSE(EnsIOUnit)
    
    call Ens_SetParams_Sat
    
    call Ens_SetParams_Seik
    
end subroutine

subroutine Ens_Finalize_Params

    integer ierror

    !deallocate(ParamArray)
    deallocate(ParamName)
    deallocate(UseParam)
    deallocate(ParamSubset)
    deallocate(SaveParam)
    
    CALL MPI_Win_free(win_ParamArray, ierror)
    
end subroutine

!     subroutine Ens_SetParams
!     
!         if (UseParam(1)) p_pu_ra(iiB1)=ParamArray(1)
!         if (UseParam(2)) p_sum(iiZ5)=ParamArray(2)
!         if (UseParam(3)) p_sum(iiZ6)=ParamArray(3)
!         if (UseParam(4)) p_qlcPPY(iiP1)=ParamArray(4)
!         if (UseParam(5)) p_qlcPPY(iiP3)=ParamArray(5)
!         if (UseParam(6)) p_qpcPPY(iiP3)=ParamArray(6)
!         if (UseParam(7)) p_qplc(iiP1)=ParamArray(7)
!         if (UseParam(8)) p_srs(iiP1)=ParamArray(8)
!         if (UseParam(9)) vsed=ParamArray(9)
!         
!     end subroutine

subroutine Ens_SetParams_trcbio(indexj, indexi)
    integer, intent(in) :: indexj, indexi

    if (UseParam(1)) p_pu_ra(iiB1)=ParamArray(1,indexj, indexi)
    if (UseParam(2)) p_sum(iiZ5)=ParamArray(2,indexj, indexi)
    if (UseParam(3)) p_sum(iiZ6)=ParamArray(3,indexj, indexi)
    if (UseParam(4)) p_qlcPPY(iiP1)=ParamArray(4,indexj, indexi)
    if (UseParam(5)) p_qlcPPY(iiP3)=ParamArray(5,indexj, indexi)
    if (UseParam(6)) p_qpcPPY(iiP3)=ParamArray(6,indexj, indexi)
    if (UseParam(7)) p_qplc(iiP1)=ParamArray(7,indexj, indexi)
    if (UseParam(8)) p_srs(iiP1)=ParamArray(8,indexj, indexi)
    
end subroutine

subroutine Ens_SetParams_trcsed(indexj, indexi)
    integer, intent(in) :: indexj, indexi
    
    if (UseParam(9)) vsed=ParamArray(9,indexj, indexi)
    
end subroutine

subroutine Ens_SetParams_Sat

    integer ierror
    
    Sat_std_log=>gl_ParamArray(11,nldj:nlej,nldi:nlei, EnsRankZero)
    Sat_std_additive=>gl_ParamArray(12,nldj:nlej,nldi:nlei, EnsRankZero)
    Sat_bias=>gl_ParamArray(13,nldj:nlej,nldi:nlei, EnsRankZero)
    win_Input_Sat=win_ParamArray
    
    CALL MPI_Win_fence(0, win_ParamArray, ierror)
    
end subroutine

subroutine Ens_SetParams_Seik

    integer ierror

    Inflation=>gl_ParamArray(10,nldj:nlej,nldi:nlei, EnsRankZero)
    win_Input_Custom=win_ParamArray
    
    CALL MPI_Win_fence(0, win_ParamArray, ierror)
    
end subroutine

    
subroutine Ens_ReadParams(prefix, opt_subset_array)
    !use Ens_Mem, &
    !    only: EnsIOUnit
    use TIME_MANAGER, &
        only: DateStart
    !use Ens_Params, &
    !    only: n_params, &
    !        ParamArray, ParamName
            
    character(LEN=*), intent(in) :: prefix
    integer, dimension(:), optional, intent(in) :: opt_subset_array
    
    integer :: indexi, ierr
    double precision :: temp
    character(LEN=200) :: filename
    logical :: existFile
    
    do indexi=1,n_params
        if (.not.(UseParam(indexi))) cycle
        if (present(opt_subset_array)) then
!             if (lwp) write(*,*) "param n. ", indexi, ", opt_subset_array: ", opt_subset_array, ", findloc: ", findloc(opt_subset_array, ParamSubset(indexi), DIM=1)
            if (findloc(opt_subset_array, ParamSubset(indexi), DIM=1)==0) cycle
        end if
        
        filename=trim(prefix)//'.'//DateStart//'.'//trim(ParamName(indexi))
        
        INQUIRE(FILE=trim(filename)//'.nc', EXIST=existFile)
        if (existFile) then
            if (lwp) write(*,*) 'Reading '//trim(filename)//'.nc'
            CALL readnc_slice_double_2d(trim(filename)//'.nc', "TRN"//trim(ParamName(indexi)), ParamArray(indexi,:,:))
            cycle
        end if
        
        INQUIRE(FILE=trim(filename)//'.txt', EXIST=existFile)
        if (existFile) then
            if (lwp) write(*,*) 'Reading '//trim(filename)//'.txt'
            open(EnsIOUnit, file=trim(filename)//'.txt', status = 'old')
                read(EnsIOUnit,*) temp
            close(EnsIOUnit)
        else
            if (lwp) write(*,*) 'Setting parameter '//trim(ParamName(indexi))//' to namelist value'
            select case (indexi)
                case (1)
                    temp=p_pu_ra(iiB1)
                case (2)
                    temp=p_sum(iiZ5)
                case (3)
                    temp=p_sum(iiZ6)
                case (4)
                    temp=p_qlcPPY(iiP1)
                case (5)
                    temp=p_qlcPPY(iiP3)
                case (6)
                    temp=p_qpcPPY(iiP3)
                case (7)
                    temp=p_qplc(iiP1)
                case (8)
                    temp=p_srs(iiP1)
                case (9)
                    temp=vsed
                case (10)
                    temp=1.0d0/sqrt(ForgettingFactor)
                case (11)
                    temp=log(1.0d0+SatMultError)
                case (12)
                    temp=SatAddError
                case (13)
                    temp=SatBias
                case default
                    if (lwp) write(*,*) "invalid parameter index in Ens_ReadParams. Aborting."
                    call mpi_barrier(glcomm,ierr)
                    call MPI_abort(glcomm, 1, ierr)
            end select
        end if
        ParamArray(indexi,:,:)=temp

    end do
    
end subroutine
    
subroutine Ens_WriteParams(DateString, prefix, opt_subset_array)

    use Ens_ParallelWriter, &
        only: PW_prepare_writing, PW_write_all
    use myalloc, &
        only: nldj, nlej, nldi, nlei
    
    character(LEN=*), intent(in) :: prefix
    character(LEN=17), intent(in) :: DateString
    integer, dimension(:), optional, intent(in) :: opt_subset_array
    
    integer :: indexi, ierr
!     character(LEN=200) :: filename
    
    do indexi=1,n_params
        if (.not.(UseParam(indexi))) cycle
        if (.not.(SaveParam(indexi))) cycle
        if (present(opt_subset_array)) then
            if (findloc(opt_subset_array, ParamSubset(indexi), DIM=1)==0) cycle
        end if
        
!         filename=trim(prefix)//'.'//DateString//'.'//trim(ParamName(indexi))//'.nc'
        
!         write(*,*) ParamName(indexi), ParamArray(indexi, nldj:nlej, nldi:nlei), "----------------------------------------------------------------------"
        call PW_prepare_writing(trim(prefix), DateString, ParamName(indexi), ParamArray(indexi, nldj:nlej, nldi:nlei), 1)
        
    end do
    
    call PW_write_all
    
end subroutine

end module
