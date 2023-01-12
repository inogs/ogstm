! developed by Simone Spada (sspada@ogs.it) at OGS

module ObsSatellite
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
        only: EnsIOUnit, &
            UseLocalObsDumping
    use Ens_Custom, &   
        only: nk_DAstate, nj_DAstate, ni_DAstate, ntra_DAstate, &
            DAMask, DAVariablesIndex
    use LocalOperations, &
        only: LocalSpace
    
    
    implicit none
    
    logical, parameter :: LocalObs=.true.
    integer, parameter :: Sat_nObs=1, Sat_n_SqrtR1=1
    double precision, parameter :: fillvalue999=-999.0d0, fillValue = 1.0d20
    
    logical :: UseSat
    double precision, dimension(:,:), allocatable :: Sat_Data, Sat_var_additive, Sat_var_log, Sat_std1_ji
    integer, dimension(4) :: SatVariablesIndex
    type(dump_container) :: satTimes
    
    type(LocalSpace) :: LocalPatch
    
    character(len=100) :: satfile_suffix, satvarname
    double precision :: SatMultError, SatAddError
    
    
contains


subroutine Sat_Namelist(filename)
            
    character(len=*) :: filename
    
    NAMELIST/Obs_Sat_setup/ UseSat, satfile_suffix, satvarname, SatMultError, SatAddError
    
    if (lwp) then
        
    end if
    
    OPEN(unit=EnsIOUnit, file=filename, status='OLD')   
        
        UseSat=.true.
        satfile_suffix='_d-OC_CNR-L3-CHL-MedOC4AD4_MULTI_1KM-MED-DT-v02.nc'
        satvarname='CHL'
        SatMultError=0.35d0
        SatAddError=0.02d0

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
            WRITE(*,*) ' '
        END IF

    CLOSE(EnsIOUnit)
    
end subroutine

Subroutine Sat_Init

    !double precision, dimension(:), POINTER, contiguous :: member_pointer
    double precision, dimension(:,:), POINTER, contiguous :: global_pointer
    
    integer :: indexi, indexj
    
    call Sat_Namelist('namelist.init')
    
    allocate(Sat_Data(jpj, jpi))
    Sat_Data=Huge(Sat_Data(1, 1))
    
    allocate(Sat_var_additive(jpj, jpi))
    Sat_var_additive=Huge(Sat_var_additive(1, 1))
    
    allocate(Sat_var_log(jpj, jpi))
    Sat_var_log=Huge(Sat_var_log(1, 1))
    
    allocate(Sat_std1_ji(nj_DAstate, ni_DAstate))
    Sat_std1_ji=Huge(Sat_std1_ji(1, 1))
    
    indexj=0
    do indexi=1, ntra_DAstate
        if (IsaSatVar(ctrcnm(DAVariablesIndex(indexi)))) then
            indexj=indexj+1
            SatVariablesIndex(indexj)=indexi
        end if
    end do
    
    satTimes%FileName = 'satTimes'
    satTimes%Name='...'
    call Load_Dump_container(satTimes)
    
    if (LocalObs) call LocalPatch%Init(1)
    
end subroutine

Subroutine Sat_Finalize
    
    deallocate(Sat_Data)
    deallocate(Sat_var_additive)
    deallocate(Sat_var_log)
    deallocate(Sat_std1_ji)
    
    call Unload_Dump_container(satTimes)
    
    if (LocalObs) call LocalPatch%Free
    
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
    call readnc_slice_double_2d(trim(VARFILE),'variance', Sat_var_additive)
    !Sat_std1_ji=1.0d0/sqrt(Sat_var_additive(nldj:nlej, nldi:nlei))
    !Sat_var_log=0.3d0**2
    Sat_var_log=log(1.0d0+SatMultError)**2
    Sat_var_additive(:,:)=SatAddError**2
    
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
                        (Sat_var_log(nldj-1+indexj,nldi-1+indexi) + log(1.0d0 + sqrt(Sat_var_additive(nldj-1+indexj,nldi-1+indexi))/tempvalue)**2) )              
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
                    Sat_std1_ji(indexj,indexi)=1.0d0/sqrt(Sat_var_log(nldj-1+indexj,nldi-1+indexi) + log(1.0d0 + sqrt(Sat_var_additive(nldj-1+indexj,nldi-1+indexi))/tempvalue)**2)                    
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
            tempvalue=Sat_Data(nldj-1+indexj,nldi-1+indexi)
            if ( (.not.(tempvalue.eq.tempvalue)).or.(tempvalue.eq.fillValue).or.(tempvalue.eq.fillvalue999).or.(DAMask(1, indexj, indexi)==0) ) then 
                Sat_Misfit(1,indexj,indexi)=0.0d0
            else
                Sat_Misfit(1, indexj,indexi)=log(Sat_Data(nldj-1+indexj, nldi-1+indexi))-ObsState(1,indexj,indexi)
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

end module
