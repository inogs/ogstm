! developed by Simone Spada (sspada@ogs.it) at OGS

module ObsSatellite_Likelihood
    use MPI
    
    use modul_param, &
        only: jpi, jpj
    use myalloc, &
        only: lwp, ctrcnm, &
            nldi, nldj, nlei, nlej
    USE ogstm_mpi_module, &
        ONLY: glcomm
    !Use IOnc, &
    !    only: readnc_slice_double_2d
    use Time_Manager, &
        only: dump_container, &
            Load_Dump_container, Unload_Dump_container
    use Ens_Mem, &
        only: UseParams, EnsRank, EnsRankZero, &
            EnsIOUnit, &
            UseLocalObsDumping, &
            Ens_shared_alloc
    use Ens_Custom, &   
        only: nk_DAstate, nj_DAstate, ni_DAstate, ntra_DAstate, &
            DAMask, DAVariablesIndex, &
            Inflation
    use LocalOperations, &
        only: LocalSpace
    
    implicit none
    
    logical, parameter :: LocalObs=.true.
    integer, parameter :: Sat_nObs=1, Sat_n_SqrtR1=1
    double precision, parameter :: fillvalue999=-999.0d0, fillValue = 1.0d20
    
    logical :: UseSat
    double precision, dimension(:,:), allocatable :: Sat_Data, Sat_std1_ji
    double precision, dimension(:,:), pointer :: Sat_std_additive, Sat_std_log, Sat_bias
    integer, dimension(4) :: SatVariablesIndex
    type(dump_container) :: satTimes
    
    type(LocalSpace) :: LocalPatch
    
    character(len=100) :: satfile_suffix, satvarname
    double precision :: SatMultError, SatAddError, SatBias
    
    !-------------------------------
    
    double precision, dimension(:,:), allocatable :: EnsCov1, TTW1T, HLTR1HL, temparray
    type(LocalSpace) :: LocalData
    double precision, dimension(:,:), allocatable :: LikelData
    double precision, dimension(:,:,:), allocatable :: ChlMean
    
    double precision, dimension(:,:), allocatable :: LikelSqrtR1
    double precision, dimension(:), allocatable :: HState
    
    double precision, pointer, contiguous, dimension(:,:,:) :: ExtendChl
    integer :: n_ExtendChl
    integer :: win_ExtendChl
    double precision, pointer, contiguous, dimension(:,:,:,:) :: gl_ExtendChl
    
    double precision :: LnDetTTW1T
    
    integer :: win_Input_Sat
    
    integer :: funcalls
    
    
contains


subroutine Sat_Namelist(filename)
            
    character(len=*) :: filename
    
    NAMELIST/Obs_Sat_setup/ UseSat, satfile_suffix, satvarname, SatMultError, SatAddError, SatBias
    
    if (lwp) then
        
    end if
    
    OPEN(unit=EnsIOUnit, file=filename, status='OLD')   
        
        UseSat=.true.
        satfile_suffix='_d-OC_CNR-L3-CHL-MedOC4AD4_MULTI_1KM-MED-DT-v02.nc'
        satvarname='CHL'
        SatMultError=0.35d0
        SatAddError=0.02d0
        SatBias=0.0d0

        REWIND(EnsIOUnit)
        READ(EnsIOUnit, Obs_Sat_setup)

        IF(lwp) THEN
            write(*,*) 'Namelist ObsSat parameters:'
            WRITE(*,*) ' '
            WRITE(*,*) ' UseSat: ', UseSat
            WRITE(*,*) ' satfile_suffix: ', trim(satfile_suffix)
            WRITE(*,*) ' satvarname: ', trim(satvarname)
            WRITE(*,*) ' SatMultError: ', SatMultError
            WRITE(*,*) ' SatAddError: ', SatAddError
            WRITE(*,*) ' SatBias: ', SatBias
            WRITE(*,*) ' '
        END IF

    CLOSE(EnsIOUnit)
    
end subroutine

Subroutine Sat_Init

    !double precision, dimension(:), POINTER, contiguous :: member_pointer
    double precision, dimension(:,:), POINTER, contiguous :: global_pointer
    
    integer :: indexi, indexj, ierror
    
    call Sat_Namelist('namelist.init')
    
    allocate(Sat_Data(jpj, jpi))
    Sat_Data=Huge(Sat_Data(1, 1))
    
    if (.not.(UseParams)) then
        allocate(Sat_std_additive(nj_DAstate, ni_DAstate))
        allocate(Sat_std_log(nj_DAstate, ni_DAstate))
        allocate(Sat_bias(nj_DAstate, ni_DAstate))
        Sat_std_additive=SatAddError
        Sat_std_log=log(1.0d0+SatMultError)
        Sat_bias=SatBias
    else
        CALL MPI_Win_fence(0, win_Input_Sat, ierror)
            if (EnsRank==EnsRankZero) then
                Sat_std_additive=SatAddError
                Sat_std_log=log(1.0d0+SatMultError)
                Sat_bias=SatBias
            end if
        CALL MPI_Win_fence(0, win_Input_Sat, ierror)
        
! write(*,*) "time: ", MPI_WTIME(), " EnsRank: ", EnsRank, ", init input_sat"
! call flush
        
    end if
    
    
    allocate(Sat_std1_ji(nj_DAstate, ni_DAstate))
    Sat_std1_ji=Huge(Sat_std1_ji(1, 1))
    
    indexj=0
    do indexi=1, ntra_DAstate
        if (IsaSatVar(trim(ctrcnm(DAVariablesIndex(indexi))))) then
            indexj=indexj+1
            SatVariablesIndex(indexj)=indexi
        end if
    end do
    
    satTimes%FileName = 'satTimes'
    satTimes%Name='...'
    call Load_Dump_container(satTimes)
    
    if (LocalObs) call LocalPatch%Init(1)
    
    call Sat_Likelihood_Init
    
end subroutine

Subroutine Sat_Finalize
    
    deallocate(Sat_Data)
    
    if (UseParams) then
        nullify(Sat_std_additive)
        nullify(Sat_std_log)
        nullify(Sat_bias)
    else
        deallocate(Sat_std_additive)
        deallocate(Sat_std_log)
        deallocate(Sat_bias)
    end if
    
    deallocate(Sat_std1_ji)
    
    call Unload_Dump_container(satTimes)
    
    if (LocalObs) call LocalPatch%Free
    
    call Sat_Likelihood_Finalize
    
end subroutine

subroutine Sat_LoadObs(DateString, nObs, n_SqrtR1)
    
    character(len=*), intent(in) :: DateString
    integer, intent(out) :: nObs, n_SqrtR1
    
    Character(len=2) :: MONTH
    Character(len=8) :: DAY
    character(LEN=1024) SATFILE, VARFILE
    integer indexi, indexj
    double precision tempvalue
    integer ierror
    
    nObs=0
    n_SqrtR1=0
    nObs=Sat_nObs
    n_SqrtR1=Sat_n_SqrtR1
    
    if ((.not.IsaDASat(DateString)).or.(.not.UseSat)) return
    
    MONTH=DateString(5:6)
    DAY  =DateString(1:8)
    
    SATFILE   = 'SATELLITE/' // DAY // trim(satfile_suffix)
    VARFILE   = 'DA_static_data/VAR_SAT/var2D.' // MONTH // '.nc'
    call readnc_slice_double_2d(trim(SATFILE),trim(satvarname), Sat_Data)
    !call readnc_slice_double_2d(trim(VARFILE),'variance', Sat_var_additive)
    !Sat_std1_ji=1.0d0/sqrt(Sat_var_additive(nldj:nlej, nldi:nlei))
    !Sat_var_log=0.3d0**2
!     Sat_std_log=log(1.0d0+SatMultError)
!     Sat_std_additive(:,:)=SatAddError

! write(*,*) "time: ", MPI_WTIME(), " EnsRank: ", EnsRank, ", before likelihood"
! call flush

    call Sat_Likelihood
    
! write(*,*) "time: ", MPI_WTIME(), " EnsRank: ", EnsRank, ", after likelihood"
! call flush
    
    if (LocalObs) then
        
        Sat_std1_ji(:,:)=1.0d0
!             do indexi=1,ni_DAstate
!                 do indexj=1,nj_DAstate
! !                     tempvalue=Sat_Data(nldj-1+indexj,nldi-1+indexi)
! !                     if ( (.not.(tempvalue.eq.tempvalue)).or.(tempvalue.eq.fillValue).or.(tempvalue.eq.fillvalue999).or.(DAMask(1, indexj, indexi)==0) ) then
!                     if (DAMask(1, indexj, indexi)==0) then
!                         Sat_std1_ji(indexj,indexi)=0.0d0
!                     else
!                         Sat_std1_ji(indexj,indexi)=1.0d0
!                     end if
!                 end do
!             end do
        
        call LocalPatch%ComputeAll(Sat_std1_ji, UseLocalObsDumping)
        
        do indexi=1,ni_DAstate
            do indexj=1,nj_DAstate
                tempvalue=Sat_Data(nldj-1+indexj,nldi-1+indexi)
                if ( (.not.(tempvalue.eq.tempvalue)).or.(tempvalue.eq.fillValue).or.(tempvalue.eq.fillvalue999).or.(DAMask(1, indexj, indexi)==0) ) then 
                    Sat_std1_ji(indexj,indexi)=0.0d0
                else
                    Sat_std1_ji(indexj,indexi)=1.0d0/sqrt( Sat_std1_ji(indexj,indexi) * &
                        (Sat_std_log(indexj,indexi)**2 + log(1.0d0 + Sat_std_additive(indexj,indexi)/tempvalue)**2) )              
                end if
            end do
        end do
        
    else
    
        if (lwp) write(*,*) 'Non-local analysis is not implemented at this moment. Aborting.'
        call MPI_abort(glcomm, 1, ierror)
    
        do indexi=1,ni_DAstate
            do indexj=1,nj_DAstate
                tempvalue=Sat_Data(nldj-1+indexj,nldi-1+indexi)
                if ( (.not.(tempvalue.eq.tempvalue)).or.(tempvalue.eq.fillValue).or.(tempvalue.eq.fillvalue999).or.(DAMask(1, indexj, indexi)==0) ) then 
                    Sat_std1_ji(indexj,indexi)=0.0d0
                else
                    Sat_std1_ji(indexj,indexi)=1.0d0/sqrt(Sat_std_log(indexj,indexi)**2 + log(1.0d0 + Sat_std_additive(indexj,indexi)/tempvalue)**2)                    
                end if
            end do
        end do
    
    end if
    
    nObs=Sat_nObs
    n_SqrtR1=Sat_n_SqrtR1
    
    if (lwp) write(*,*) 'Loaded '//trim(SATFILE)
    
end subroutine

Function Sat_Misfit(ObsState)
    double precision, dimension(Sat_nObs,nj_DAstate, ni_DAstate), intent(in) :: ObsState
    double precision, dimension(Sat_nObs,nj_DAstate, ni_DAstate) :: Sat_Misfit
    
    integer indexi, indexj
    double precision tempvalue
    
    if (.not.UseSat) return
    
    !Sat_Misfit(1,:,:)=Sat_Data(nldj:nlej, nldi:nlei)-ObsState(1,:,:)
    do indexi=1,ni_DAstate
        do indexj=1,nj_DAstate
            tempvalue=Sat_Data(nldj-1+indexj, nldi-1+indexi)
            if ( (.not.(tempvalue.eq.tempvalue)).or.(tempvalue.eq.fillValue).or.(tempvalue.eq.fillvalue999).or.(DAMask(1, indexj, indexi)==0) ) then 
                Sat_Misfit(1,indexj,indexi)=0.0d0
            else
                Sat_Misfit(1, indexj,indexi)=log(tempvalue)-ObsState(1,indexj,indexi)-Sat_bias(indexj,indexi)
            end if
        end do
    end do
    
end function

function Sat_H(State)
    
    double precision, dimension(nk_DAstate,nj_DAstate,ni_DAstate,ntra_DAstate), intent(in) :: State
    double precision, dimension(Sat_nObs, nj_DAstate,ni_DAstate) :: Sat_H
    
    integer indexi, indexj
    double precision tempvalue
    
    if (.not.UseSat) return
    
    !Sat_H=log(sum(exp(State(1,:,:,SatVariablesIndex), dim=3)))
    
    do indexi=1,ni_DAstate
        do indexj=1,nj_DAstate
            tempvalue=Sat_Data(nldj-1+indexj,nldi-1+indexi)
            if ( (.not.(tempvalue.eq.tempvalue)).or.(tempvalue.eq.fillValue).or.(tempvalue.eq.fillvalue999).or.(DAMask(1, indexj, indexi)==0) ) then 
                Sat_H(1,indexj,indexi)=0.0d0
            else                    
                Sat_H(1,indexj,indexi)=log(sum(exp(State(1,indexj,indexi,SatVariablesIndex))))              
            end if
        end do
    end do

end function

function Sat_SqrtR1(HLi)
    double precision, dimension(Sat_nObs, nj_DAstate, ni_DAstate), intent(in) :: HLi
    double precision, dimension(Sat_n_SqrtR1, nj_DAstate, ni_DAstate) :: Sat_SqrtR1
    
    if (.not.UseSat) return
    
    Sat_SqrtR1(1, :,:)=Sat_std1_ji*HLi(1,:,:)
    
end function

function IsaSatVar(name)
    implicit none
    
    character(LEN=*), intent(in) :: name
    logical :: IsaSatVar
    
    if ((name.eq."P1l").or.(name.eq."P2l").or.(name.eq."P3l").or.(name.eq."P4l")) then
        IsaSatVar=.true.
    else
        IsaSatVar=.false.
    end if
        
end function

LOGICAL FUNCTION IsaDASat(datestring)
            
    CHARACTER(LEN=17), INTENT(IN) :: datestring
    ! LOCAL
    INTEGER I

    IsaDASat = .false.
    
    DO I=1, satTimes%N
        if (datestring.eq.satTimes%TimeStrings(I)) then
            IsaDASat = .true.
            return
        endif
    ENDDO
    
end function

Subroutine Sat_Likelihood_Init

    use myalloc, &
        only: lwp, myrank, mysize, mycomm
    use Ens_Mem, &
        only: EnsRankZero, myrankZero, &
            EnsDebug, EnsRank, EnsSize, &
            LocalRange, EnsWeights
    use Ens_Custom, &   
        only: nk_DAstate, nj_DAstate, ni_DAstate, ntra_DAstate, &
            DAstate_kjit, win_DAstate, n_DAstate, gl_DAstate, gl_DAstate_kjitn, &
            DAMask
    use Ens_Utilities, &
        only: Ens_ReduceMeanAndBase
    use Algebra, &
        only: Init_Algebra, Finalize_Algebra, SymChangeBase
    use LBFGSB, &
        only: LBFGSB_Init, LBFGSB_Finalize, LBFGSB_OPT
        
    double precision, dimension(:), POINTER, contiguous :: member_pointer
    double precision, dimension(:,:), POINTER, contiguous :: global_pointer
    
    integer indexi
    
    n_ExtendChl=(nj_DAstate+2*LocalRange)*(ni_DAstate+2*LocalRange)*size(SatVariablesIndex)
    call Ens_shared_alloc(n_ExtendChl, member_pointer, global_pointer, win_ExtendChl)
    ExtendChl(1:size(SatVariablesIndex), 1-LocalRange:nj_DAstate+LocalRange, 1-LocalRange:ni_DAstate+LocalRange)=>member_pointer
    gl_ExtendChl(1:size(SatVariablesIndex), 1-LocalRange:nj_DAstate+LocalRange, 1-LocalRange:ni_DAstate+LocalRange, 0:EnsSize-1)=>global_pointer
    ExtendChl  = huge(ExtendChl(1,1,1))

    allocate(TTW1T(EnsSize-1,EnsSize-1))
    do indexi=1, EnsRankZero
        TTW1T(:,indexi)=-1.0d0
        TTW1T(indexi,indexi)=1.0d0/EnsWeights(indexi-1)-1.0d0
    end do
    do indexi=EnsRankZero+1, EnsSize-1
        TTW1T(:,indexi)=-1.0d0
        TTW1T(indexi,indexi)=1.0d0/EnsWeights(indexi)-1.0d0
    end do
    LnDetTTW1T=sum(log(EnsWeights))-log(EnsWeights(EnsRankZero))+log(1-sum(EnsWeights) + EnsWeights(EnsRankZero))
    
    allocate(EnsCov1(EnsSize-1,EnsSize-1))
    EnsCov1=Huge(EnsCov1(1,1))
    
    allocate(HLTR1HL(EnsSize-1,0:EnsSize-1))
    HLTR1HL=Huge(HLTR1HL(1, 0))
    
    call LocalData%Init(1)
    
    allocate(LikelData(nj_DAstate, ni_DAstate))
    LikelData=Huge(LikelData(1,1))
    
    allocate(ChlMean(1:size(SatVariablesIndex), 1-LocalRange:nj_DAstate+LocalRange, 1-LocalRange:ni_DAstate+LocalRange))
    ChlMean=Huge(ChlMean(1,1,1))
    
    allocate(LikelSqrtR1(1-LocalRange:1+LocalRange, 1-LocalRange:1+LocalRange))
    LikelSqrtR1=Huge(LikelSqrtR1(1,1))
    
    allocate(HState(0:EnsSize-1))
    HState=Huge(HState(0))
    
    call LBFGSB_Init

end subroutine

subroutine Sat_likelihood_Finalize
    
    use myalloc, &
        only: lwp, myrank, mysize, mycomm
    use Ens_Mem, &
        only: EnsRankZero, myrankZero, &
            EnsDebug, EnsRank, EnsSize, &
            LocalRange, EnsWeights
    use Ens_Custom, &   
        only: nk_DAstate, nj_DAstate, ni_DAstate, ntra_DAstate, &
            DAstate_kjit, win_DAstate, n_DAstate, gl_DAstate, gl_DAstate_kjitn, &
            DAMask
    use Ens_Utilities, &
        only: Ens_ReduceMeanAndBase
    use Algebra, &
        only: Init_Algebra, Finalize_Algebra, SymChangeBase
    use LBFGSB, &
        only: LBFGSB_Init, LBFGSB_Finalize, LBFGSB_OPT

    integer ierror

    deallocate(TTW1T)
    deallocate(EnsCov1)
    deallocate(LikelData)
    deallocate(LikelSqrtR1)
    deallocate(ChlMean)
    deallocate(HState)
    deallocate(HLTR1HL)
    
    call LocalData%Free
    
    CALL MPI_Win_free(win_ExtendChl, ierror)
    
    call LBFGSB_Finalize

end subroutine

subroutine Sat_Likelihood

    use myalloc, &
        only: lwp, myrank, mysize, mycomm
    use Ens_Mem, &
        only: EnsRankZero, myrankZero, &
            EnsDebug, EnsRank, EnsSize, &
            LocalRange, EnsWeights
    use Ens_Custom, &   
        only: nk_DAstate, nj_DAstate, ni_DAstate, ntra_DAstate, &
            DAstate_kjit, win_DAstate, n_DAstate, gl_DAstate, gl_DAstate_kjitn, &
            DAMask
    use Ens_Utilities, &
        only: Ens_ReduceMeanAndBase
    use Algebra, &
        only: Init_Algebra, Finalize_Algebra, SymChangeBase
    use LBFGSB, &
        only: LBFGSB_Init, LBFGSB_Finalize, LBFGSB_OPT
    
    integer :: indexi, indexj, c
    integer :: ierror
    double precision :: tempvalue
    double precision, dimension(3) :: input_parameters, lower_bound, upper_bound
    integer, dimension(2) :: other_parameters
    
    lower_bound(1)   = 1.0d-3
    lower_bound(2)   = 1.0d-3
    lower_bound(3)   = 0.0d0
    
    upper_bound(1)   = 3.0d0
    upper_bound(2)   = 5.0d0
    upper_bound(3)   = 5.0d0
    
! write(*,*) "time: ", MPI_WTIME(), " EnsRank: ", EnsRank, ", myrank: ", myrank, ", start likelihood"
! call flush
    
    do indexi=1, size(SatVariablesIndex)
    
!         call LocalData%SendRecive(DAstate_kjit(1,:,:,SatVariablesIndex(indexi)))
        allocate(temparray(nj_DAstate, ni_DAstate))
        temparray(:,:)=DAstate_kjit(1,:,:,SatVariablesIndex(indexi))
        call LocalData%SendRecive(temparray)
        deallocate(temparray)
        
        ExtendChl(indexi, :, :) = LocalData%LocalPatch(1,:,:)
    end do
    
! write(*,*) "time: ", MPI_WTIME(), " EnsRank: ", EnsRank, ", myrank: ", myrank, ", before likeldata"
! call flush
    
    LikelData(:, :) = DAMask(1, :, :)
    
! write(*,*) "time: ", MPI_WTIME(), " EnsRank: ", EnsRank, ", myrank: ", myrank, ", after likeldata"
! call flush
    
    do indexi=1,ni_DAstate
    
! write(*,*) "time: ", MPI_WTIME(), " EnsRank: ", EnsRank, ", myrank: ", myrank, ", indexi: ", indexi

        do indexj=1,nj_DAstate

! write(*,*) "time: ", MPI_WTIME(), " EnsRank: ", EnsRank, ", myrank: ", myrank, ", i: ", indexi, ", j: ", indexj, ", before tempvalue"

            tempvalue=Sat_Data(nldj-1+indexj,nldi-1+indexi)

! write(*,*) "time: ", MPI_WTIME(), " EnsRank: ", EnsRank, ", myrank: ", myrank, ", i: ", indexi, ", j: ", indexj, ", tempvalue: ", tempvalue 

            if ( (.not.(tempvalue.eq.tempvalue)).or.(tempvalue.eq.fillValue).or.(tempvalue.eq.fillvalue999).or.(LikelData(indexj, indexi)<0.5d0).or.(tempvalue<=0.0d0) ) then 
!             if ( (.not.(tempvalue.eq.tempvalue)).or.(tempvalue.eq.fillValue).or.(tempvalue.eq.fillvalue999).or.(DAMask(1, indexj, indexi)==0).or.(tempvalue<=0.0d0) ) then 
                LikelData(indexj,indexi)=fillvalue999
            else
            
! write(*,*) "time: ", MPI_WTIME(), " EnsRank: ", EnsRank, ", myrank: ", myrank, ", i: ", indexi, ", j: ", indexj, ", before satdata"            
            
                LikelData(indexj,indexi)=tempvalue
                
! write(*,*) "time: ", MPI_WTIME(), " EnsRank: ", EnsRank, ", myrank: ", myrank, ", i: ", indexi, ", j: ", indexj, ", after satdata"

            end if
        end do
    end do
    
! write(*,*) "time: ", MPI_WTIME(), " EnsRank: ", EnsRank, ", myrank: ", myrank, ", before likeldata sendrecive"
! call flush
    
    call LocalData%SendRecive(LikelData, opt_FillValue=fillvalue999)
    
! write(*,*) "time: ", MPI_WTIME(), " EnsRank: ", EnsRank, ", myrank: ", myrank, ", before win_ExtendChl"
! call flush
    
    CALL MPI_Win_fence(0, win_ExtendChl, ierror)
        
        ChlMean=0.0d0
        do indexi=0,EnsSize-1
            ChlMean=ChlMean+gl_ExtendChl(:,:,:,indexi)*EnsWeights(indexi)
        end do
        
        ExtendChl=ExtendChl - ChlMean
        
! write(*,*) "time: ", MPI_WTIME(), " EnsRank: ", EnsRank, ", myrank: ", myrank, ", before win_Input_Sat"
! call flush
        
    CALL MPI_Win_fence(0, win_ExtendChl, ierror)
    CALL MPI_Win_fence(0, win_Input_Sat, ierror)
    
! write(*,*) "time: ", MPI_WTIME(), " EnsRank: ", EnsRank, ", myrank: ", myrank, ", after win_Input_Sat"
! call flush
        
        c=-1
        do indexi=1,ni_DAstate
            do indexj=1,nj_DAstate
                tempvalue=LocalData%LocalPatch(1,indexj,indexi)
                if (tempvalue.eq.fillvalue999) then
                    if (EnsRank==EnsRankZero) Inflation(indexj,indexi)=1.0d0
                    cycle
                end if
                c=c+1
                if (EnsRank/=mod(c,EnsSize)) cycle
                
                input_parameters(1)=Sat_std_log(indexj,indexi)*10.0d0
                input_parameters(2)=Sat_std_additive(indexj,indexi)*100.0d0
                input_parameters(3)=(Inflation(indexj,indexi)-1.0d0)*100.0d0
                
                other_parameters(1)=indexi
                other_parameters(2)=indexj
                
                funcalls=0
                
! write(*,*) "EnsRank: ", EnsRank, ", myrank: ", myrank, ", before LBFGSB_OPT"
! call flush
                
                call LBFGSB_OPT(MinFun, other_parameters, lower_bound, upper_bound, input_parameters, tempvalue, ierror)
                
! write(*,*) "EnsRank: ", EnsRank, ", myrank: ", myrank, ", after LBFGSB_OPT"
! call flush

                
                if (ierror==0) then        
!                     write(*,*) "myrank: ", myrank, ", i: ",indexi, ", indexj: ",indexj, ", funcalls: ", funcalls, ", f: ", tempvalue, ", inputold: ", Sat_std_log(nldj-1+indexj,nldi-1+indexi), Sat_std_additive(nldj-1+indexj,nldi-1+indexi), Inflation(nldj-1+indexj,nldi-1+indexi), ", inputnew: ", input_parameters(1)*1.0d-1, input_parameters(2)*1.0d-2, 1.0d0+input_parameters(3)*1.0d-2  
                    
                    Sat_std_log(indexj,indexi)=input_parameters(1)*1.0d-1
                    Sat_std_additive(indexj,indexi)=input_parameters(2)*1.0d-2
                    Inflation(indexj,indexi)=1.0d0+input_parameters(3)*1.0d-2            
                else
!                     write(*,*) "myrank: ", myrank, ", i: ",indexi, ", indexj: ",indexj, ", funcalls: ", funcalls, ", ierror: ", ierror
                end if
            end do
        end do
        
    CALL MPI_Win_fence(0, win_ExtendChl, ierror)
    CALL MPI_Win_fence(0, win_Input_Sat, ierror)
    
end subroutine
    
subroutine Sat_Likelihood_Score(indexj, indexi, input_parameters, Score, ierror)
    use myalloc, &
        only: lwp, myrank, mysize, mycomm
    use Ens_Mem, &
        only: EnsRankZero, myrankZero, &
            EnsDebug, EnsRank, EnsSize, &
            LocalRange, EnsWeights
    use Ens_Custom, &   
        only: nk_DAstate, nj_DAstate, ni_DAstate, ntra_DAstate, &
            DAstate_kjit, win_DAstate, n_DAstate, gl_DAstate, gl_DAstate_kjitn, &
            DAMask
    use Ens_Utilities, &
        only: Ens_ReduceMeanAndBase
    use Algebra, &
        only: Init_Algebra, Finalize_Algebra, SymChangeBase
    use LBFGSB, &
        only: LBFGSB_Init, LBFGSB_Finalize, LBFGSB_OPT
            
    integer, intent(in) :: indexi, indexj
    double precision, dimension(3), intent(in) :: input_parameters 
    double precision, intent(out) :: Score
    integer, intent(out) :: ierror
            
    double precision, parameter :: pi=3.1415926535d0
    integer ::  indexk, indexw, indexd, itemp
    double precision :: tempvalue, LnDetSqrtR1, SqrtR1ij
    double precision :: LnMultErr, AddErr, Infl, Bias
    
    LnDetSqrtR1=0.0d0
    HLTR1HL=0.0d0
    Score=0.0d0
    
    LnMultErr=input_parameters(1)*1.0d-1
    AddErr=input_parameters(2)*1.0d-2
    Infl=input_parameters(3)*1.0d-2 + 1.0d0
    !Bias=input_parameters(4) !Sat_bias(nldj-1+indexj,nldi-1+indexi)
    Bias=0.0d0
    
    do indexw=-LocalRange,LocalRange
        itemp=floor(sqrt((0.5d0+LocalRange)**2-indexw**2))
        do indexk=-itemp,itemp
            tempvalue=LocalData%LocalPatch(1,indexj+indexk,indexi+indexw)
            if (tempvalue.eq.fillvalue999) cycle
            
            !SqrtR1ij=1.0d0/sqrt( (LocalRange+0.5d0)**2*pi * &
            !    (Sat_var_log(nldj-1+indexj,nldi-1+indexi) + log(1.0d0 + sqrt(Sat_var_additive(nldj-1+indexj,nldi-1+indexi))/tempvalue)**2) )
            
!             write(*,*) "myrank: ", myrank, ", i: ",indexi, ", indexj: ",indexj, ", LnMultErr: ", LnMultErr, ", AddErr: ",AddErr, ", Infl: ",Infl, ", Bias: ", Bias, ", tempvalue: ", tempvalue
                
            SqrtR1ij=1.0d0/sqrt( (LocalRange+0.5d0)**2*pi * (LnMultErr**2 + log(1.0d0 + AddErr/tempvalue)**2) )
            
!             write(*,*) "myrank: ", myrank, ", i: ",indexi, ", indexj: ",indexj, ", SqrtR1ij: ", SqrtR1ij
            
            LnDetSqrtR1=LnDetSqrtR1+log(SqrtR1ij)
                
            do indexd=0, EnsSize-1
                HState(indexd) = log(sum(exp(gl_ExtendChl(:,indexj+indexk,indexi+indexw, indexd)*Infl + &
                    ChlMean(:,indexj+indexk,indexi+indexw) )))
            end do
            
!             write(*,*) "myrank: ", myrank, ", i: ",indexi, ", indexj: ",indexj, ", HState1: ", HState
            
            HState(EnsRankZero)= DOT_PRODUCT(HState, EnsWeights)
            
            do indexd=0, EnsSize-1
                if (indexd==EnsRankZero) cycle
                HState(indexd) = HState(indexd)-HState(EnsRankZero)
            end do
            
!             write(*,*) "myrank: ", myrank, ", i: ",indexi, ", indexj: ",indexj, ", HState2: ", HState
            
            HState(EnsRankZero)=log(tempvalue) - HState(EnsRankZero) - Bias
            
            HState(:) = HState*SqrtR1ij
            
!             write(*,*) "myrank: ", myrank, ", i: ",indexi, ", indexj: ",indexj, ", HState3: ", HState
            
            do indexd=0, EnsSize-1
                HLTR1HL(1:EnsRankZero, indexd) = HLTR1HL(1:EnsRankZero, indexd) + HState(0:EnsRankZero-1)*HState(indexd)
                HLTR1HL(EnsRankZero+1:EnsSize-1, indexd) = HLTR1HL(EnsRankZero+1:EnsSize-1, indexd) + HState(EnsRankZero+1:EnsSize-1)*HState(indexd)
            end do
            
!             write(*,*) "myrank: ", myrank, ", i: ",indexi, ", indexj: ",indexj, ", HLTR1HL: ", HLTR1HL
            
            Score = Score+HState(EnsRankZero)**2
            
        end do
    end do
    
    EnsCov1(:,1:EnsRankZero)=TTW1T(:,1:EnsRankZero)+HLTR1HL(:, 0:EnsRankZero-1)
    EnsCov1(:,EnsRankZero+1:EnsSize-1)=TTW1T(:,EnsRankZero+1:EnsSize-1)+HLTR1HL(:, EnsRankZero+1:EnsSize-1)
    
    call SymChangeBase(EnsCov1,ierror, LnDet_opt=tempvalue)
    if (ierror/=0) then
        write(*,*) "Rank: ", myrank, ", indexi: ", indexi, ", indexj: ", indexj, " -> SymChangeBase failed in Sat_Likelihood_Score."
        return
    end if
    
    HState(1:EnsSize-1)=matmul(EnsCov1,HLTR1HL(:,EnsRankZero))
    
!     write(*,*) "myrank: ", myrank, ", i: ",indexi, ", indexj: ",indexj, ", Score after loop: ", Score, ", dot: ", DOT_PRODUCT(HState(1:EnsSize-1), HState(1:EnsSize-1)), ", tempvalue: ", tempvalue, ", LnDetSqrtR1: ", LnDetSqrtR1 
    
    Score = Score - DOT_PRODUCT(HState(1:EnsSize-1), HState(1:EnsSize-1)) + tempvalue - LnDetTTW1T - 2*LnDetSqrtR1
    
!     write(*,*) "myrank: ", myrank, ", i: ",indexi, ", indexj: ",indexj, ", Score final: ", Score
!     write(*,*) "--------------------------------------------------------------------------------------"
    
end subroutine

subroutine MinFun(input_parameters, other_parameters, f, df, ierror)
    use myalloc, &
        only: lwp, myrank, mysize, mycomm
    use Ens_Mem, &
        only: EnsRankZero, myrankZero, &
            EnsDebug, EnsRank, EnsSize, &
            LocalRange, EnsWeights
    use Ens_Custom, &   
        only: nk_DAstate, nj_DAstate, ni_DAstate, ntra_DAstate, &
            DAstate_kjit, win_DAstate, n_DAstate, gl_DAstate, gl_DAstate_kjitn, &
            DAMask
    use Ens_Utilities, &
        only: Ens_ReduceMeanAndBase
    use Algebra, &
        only: Init_Algebra, Finalize_Algebra, SymChangeBase
    use LBFGSB, &
        only: LBFGSB_Init, LBFGSB_Finalize, LBFGSB_OPT

    integer, parameter :: n=3

    double precision, dimension(n), intent(in) :: input_parameters
    integer, dimension(2), intent(in) :: other_parameters
    
    double precision, intent(out) :: f
    double precision, dimension(n), intent(out) :: df
    integer, intent(out) :: ierror
    
    double precision, parameter :: delta=1.0d-2
    
    double precision :: f2
    double precision, dimension(3) :: new_param
    
    integer indexi, indexj, indexk
    
    funcalls=funcalls+1
    
    !LnMultErr=input_parameters(1)
    !AddErr=input_parameters(2)
    !Infl=input_parameters(3)
    !Bias=input_parameters(4) !Sat_bias(nldj-1+indexj,nldi-1+indexi)
    
    indexi=other_parameters(1)
    indexj=other_parameters(2)
    
!test
! f=sum((input_parameters-1.0d0)**2)
! ierror=0
!!!!!!!!!!!!!
    call Sat_Likelihood_Score(indexj, indexi, input_parameters, f, ierror)
    if (ierror/=0) return
    
    do indexk=1,n
        new_param=input_parameters
        new_param(indexk)=input_parameters(indexk)+delta
!test
! f2=sum((new_param-1.0d0)**2)
! ierror=0
!!!!!!!!!!!!!
        call Sat_Likelihood_Score(indexj, indexi, new_param, f2, ierror)
        if (ierror/=0) return
        
        df(indexk)=(f2-f)/delta
    end do
!     write(*,*) "-----------------------------------------------------------------------------------------------------"
!     write(*,*) "myrank: ", myrank, ", i: ",indexi, ", indexj: ",indexj, ", f: ", f, ", df: ", df, ", funcalls: ", funcalls, ", input_parameters: ", input_parameters
!     write(*,*) "-----------------------------------------------------------------------------------------------------"
    
end subroutine



end module


