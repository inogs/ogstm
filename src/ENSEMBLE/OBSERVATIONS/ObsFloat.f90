! developed by Simone Spada (sspada@ogs.it) at OGS

module ObsFloat
    use MPI
    use netcdf
    
    use modul_param, &
        only: jpi, jpj
    use myalloc, &
        only: lwp, ctrcnm, &
            nldi, nldj, nlei, nlej, &
            gphit, glamt, gdept, mbathy
    USE ogstm_mpi_module, &
        ONLY: glcomm
    !Use IOnc, &
    !    only: readnc_slice_double_2d
    use Time_Manager, &
        only: dump_container, &
            Load_Dump_container, Unload_Dump_container
    use Ens_Mem, &
        only: EnsIOUnit
    use Ens_Custom, &   
        only: nk_DAstate, nj_DAstate, ni_DAstate, ntra_DAstate, &
            DAMask, DAVariablesIndex    
    
    implicit none
    
    logical :: UseFloat
    integer :: Float_nObs, Float_n_SqrtR1
    integer, dimension(4) :: FloatVariablesIndex
    type(dump_container) :: FloatTimes
    
    type :: FloatType
        integer nLevels
        character(len=1024) :: FloatFILE, VarName
        character(len=17) :: FloatDate
        double precision :: FloatLat, FloatLon
        integer :: i, j
        real, dimension(:), allocatable :: tracer, depth, qc, std1, k
        !integer, dimension(:), allocatable :: k
    contains
        procedure Init => FloatType_Init
        procedure Free => FloatType_Free
    end type
    
    integer :: nFloats
    type(FloatType), dimension(:), allocatable :: Floats
    
    double precision :: FloatMultError, FloatAddError
    
    
contains

subroutine Float_Namelist(filename)
            
    character(len=*) :: filename
    
    NAMELIST/Obs_Float_setup/ UseFloat, FloatMultError, FloatAddError
    
    if (lwp) then
        
    end if
    
    OPEN(unit=EnsIOUnit, file=filename, status='OLD')
        
        UseFloat=.true.
        FloatMultError=0.35d0
        FloatAddError=0.02d0

        REWIND(EnsIOUnit)
        READ(EnsIOUnit, Obs_Float_setup)

        IF(lwp) THEN
            write(*,*) 'Namelist ObsFloat parameters:'
            WRITE(*,*) ' '
            WRITE(*,*) ' UseFloat: ', UseFloat
            WRITE(*,*) ' FloatMultError: ', FloatMultError
            WRITE(*,*) ' FloatAddError: ', FloatAddError
            WRITE(*,*) ' '
        END IF

    CLOSE(EnsIOUnit)
    
end subroutine

subroutine FloatType_Init(this, month)
    class(FloatType), intent(inout) :: this
    Character(len=2), intent(in) :: month
    
    integer ncid, stat, objid
    real, dimension(12) :: errorfloat
    
    errorfloat = (/0.0690, 0.0969, 0.0997, 0.0826, 0.0660, 0.0500, 0.0360, 0.0140, 0.0320, 0.0390, 0.0340, 0.0490/)
    
    stat = nf90_open(trim(this%FloatFILE), nf90_nowrite, ncid)  
    call nf90error(trim(this%FloatFILE), stat)
    
    stat = nf90_inq_dimid(ncid, 'n'//trim(this%VarName), objid)
    call nf90error(trim(this%FloatFILE), stat)
    
    stat=nf90_inquire_dimension(ncid, objid, len=this%nLevels)
    call nf90error(trim(this%FloatFILE), stat)
    
    call getvar(ncid, this%nLevels, trim(this%VarName), trim(this%FloatFILE), this%tracer)
    call getvar(ncid, this%nLevels, 'PRES_'//trim(this%VarName), trim(this%FloatFILE), this%depth)
    call getvar(ncid, this%nLevels, trim(this%VarName)//'_QC', trim(this%FloatFILE), this%qc)
    
    stat = nf90_close(ncid)       
    call nf90error(trim(this%FloatFILE), stat)
    
    allocate(this%std1(this%nLevels))
    
    !read(month,*) ncid
    !this%std=errorfloat(ncid)
    
    this%std1 = 1.0d0 / sqrt( this%nLevels * (log(1.0d0+FloatMultError)**2 +  log(1.0d0 + FloatAddError/this%tracer)**2) )    
    
end subroutine

subroutine nf90error(filename, stat)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: stat
    
    integer ierror
    
    if (stat==nf90_noerr) return
    
    write(*,*) 'error with file:', filename
    write(*,*) nf90_strerror(stat)
    
    call MPI_abort(glcomm, 1, ierror)
    
end subroutine

subroutine getvar(ncid, nLevels, varname, filename, var)
    integer, intent(in) :: ncid, nLevels
    character(len=*), intent(in) :: varname, filename
    real, dimension(:), allocatable, intent(out) :: var
    
    integer :: objid, stat
    
    allocate(var(nLevels))
    
    stat = nf90_inq_varid (ncid, trim(varname), objid)
    call nf90error(trim(filename), stat)
    
    stat = nf90_get_var (ncid,objid,var)
    call nf90error(trim(filename), stat)
    
end subroutine

subroutine FloatType_Free(this)
    class(FloatType), intent(inout) :: this
    
    if (allocated(this%tracer)) deallocate(this%tracer)
    if (allocated(this%depth)) deallocate(this%depth)
    if (allocated(this%qc)) deallocate(this%qc)
    if (allocated(this%std1)) deallocate(this%std1)
    if (allocated(this%k)) deallocate(this%k)

end subroutine

Subroutine Float_Init

    !double precision, dimension(:), POINTER, contiguous :: member_pointer
    double precision, dimension(:,:), POINTER, contiguous :: global_pointer
    
    integer :: indexi, indexj
    
    call Float_Namelist('namelist.init')
    
    indexj=0
    do indexi=1, ntra_DAstate
        if (IsaFloatVar(ctrcnm(DAVariablesIndex(indexi)))) then
            indexj=indexj+1
            FloatVariablesIndex(indexj)=indexi
        end if
    end do
    
    FloatTimes%FileName = 'floatTimes'
    FloatTimes%Name='...'
    call Load_Dump_container(floatTimes)
    
    nfloats=0
    
end subroutine

Subroutine Float_Finalize

    integer indexi

    if (allocated(Floats)) then
        do indexi=1,nFloats
            call Floats(indexi)%Free
        end do
        deallocate(Floats)
    end if
    
    call Unload_Dump_container(FloatTimes)
    
end subroutine

subroutine Float_LoadObs(DateString, nObs, n_SqrtR1)
    
    character(len=*), intent(in) :: DateString
    integer, intent(out) :: nObs, n_SqrtR1
    
    Character(len=2) :: MONTH
    Character(len=8) :: DAY
    character(LEN=1024) VARFILE, FloatIndexFile
    character(len=4) :: varname
    integer indexi, indexj, indexk
    double precision tempvalue
    integer ierror
    
    character(len=1024) :: line, FloatFILE
    character(len=17) :: FloatDate
    double precision :: FloatLat, FloatLon
    double precision :: dLat, dLon
    
    nObs=0
    n_SqrtR1=0
    Float_nObs=nObs
    Float_n_SqrtR1=n_SqrtR1
    
    if ((.not.IsaDAFloat(DateString)).or.(.not.UseFloat)) return
    
    if (allocated(Floats)) then
        do indexi=1,nFloats
            call Floats(indexi)%Free
        end do
        deallocate(Floats)
    end if
    
    MONTH=DateString(5:6)
    DAY  =DateString(1:8)
    
    FloatIndexFile   = 'FLOATS/Float_Index.txt'   
    varname='CHLA'
    
    !VARFILE   = 'DA_static_data/VAR_Float/var2D.' // MONTH // '.nc'
    
    dLat=(gphit(2,1)-gphit(1,1))*0.5d0
    dLon=(glamt(2,1)-glamt(1,1))*0.5d0
        
    OPEN(UNIT=EnsIOUnit,file=FloatIndexFile, status='old', position='rewind')
        do indexj=1,2
            nFloats=0
            DO 
                read(EnsIOUnit, '(a)', iostat=ierror) line
                if (ierror/=0) exit
                
                if (Index(line, varname)==0) cycle
                
                indexi=index(line,',')
                FloatFILE=line(1:indexi-1)
                line=trim(line(indexi+1:len(line)))
                
                indexi=index(line,',')
                indexi=index(line(indexi+1:len(line)),',')+indexi
                read(line(1:indexi-1),*) FloatLat, FloatLon
                line=line(indexi+1:len(line))
                
                indexi=index(line,',')
                FloatDate=line(1:indexi-1)
                line=trim(line(indexi+1:len(line)))
                
                if (FloatDate(1:8)/=DAY) cycle
                
                nFloats = nFloats+1
                
                if (indexj==1) cycle
                
                Floats(nFloats)%FloatFILE='FLOATS/'//trim(FloatFILE)
                Floats(nFloats)%FloatDate=FloatDate
                Floats(nFloats)%VarName=trim(varname)
                Floats(nFloats)%FloatLat=FloatLat
                Floats(nFloats)%FloatLon=FloatLon
                
                call Floats(nFloats)%Init(MONTH)
                
                nObs=nObs + Floats(nFloats)%nLevels
                n_SqrtR1=n_SqrtR1 + Floats(nFloats)%nLevels
                
                Floats(nFloats)%i=0
                Floats(nFloats)%j=0
                if (FloatLat<gphit(nldj,1)-dLat .or. gphit(nlej,1)+dLat<FloatLat) cycle
                if (FloatLon<glamt(1,nldi)-dLon .or. glamt(1,nlei)+dLon<FloatLon) cycle
                Floats(nFloats)%j=minloc(abs(FloatLat-gphit(nldj:nlej,1)), dim=1)
                Floats(nFloats)%i=minloc(abs(FloatLon-glamt(1,nldi:nlei)), dim=1)
                
                allocate(Floats(nFloats)%k(Floats(nFloats)%nLevels))
                indexk=1
                indexi=1
                
#ifdef gdept1d
                do while (Floats(nFloats)%depth(indexi)<gdept(1))
                    indexi=indexi+1
                end do
                Floats(nFloats)%k(1:indexi-1)=1.0d0
                do while (indexi <= Floats(nFloats)%nLevels)
                    if (Floats(nFloats)%depth(indexi)<gdept(indexk+1)) then
                        tempvalue=(Floats(nFloats)%depth(indexi)-gdept(indexk)) / (gdept(indexk+1)-gdept(indexk))
                        Floats(nFloats)%k(indexi)=indexk+tempvalue
                        indexi=indexi+1
                    elseif (indexk==mbathy(Floats(nFloats)%j,Floats(nFloats)%i)) then
                        Floats(nFloats)%k(indexi:Floats(nFloats)%nLevels)=indexk
                        exit
                    else
                        indexk=indexk+1
                    end if
                end do
#else
                do while (Floats(nFloats)%depth(indexi)<gdept(1,Floats(nFloats)%j,Floats(nFloats)%i))
                    indexi=indexi+1
                end do
                Floats(nFloats)%k(1:indexi-1)=1.0d0
                do while (indexi <= Floats(nFloats)%nLevels)
                    if (Floats(nFloats)%depth(indexi)<gdept(indexk+1,Floats(nFloats)%j,Floats(nFloats)%i)) then
                        tempvalue=(Floats(nFloats)%depth(indexi)-gdept(indexk,Floats(nFloats)%j,Floats(nFloats)%i)) / &
                            (gdept(indexk+1,Floats(nFloats)%j,Floats(nFloats)%i)-gdept(indexk,Floats(nFloats)%j,Floats(nFloats)%i))
                        Floats(nFloats)%k(indexi)=indexk+tempvalue
                        indexi=indexi+1
                    elseif (indexk==mbathy(Floats(nFloats)%j,Floats(nFloats)%i)) then
                        Floats(nFloats)%k(indexi:Floats(nFloats)%nLevels)=indexk
                        exit
                    else
                        indexk=indexk+1
                    end if
                end do
#endif
                
            ENDDO
            
            if (nFloats==0) exit
            
            if (indexj==1) then
                REWIND(EnsIOUnit)
                allocate(Floats(nFloats))
            end if
        end do
        
    CLOSE(EnsIOUnit)
        
    Float_nObs=nObs
    Float_n_SqrtR1=n_SqrtR1
    
    if (lwp) write(*,*) 'Loaded ', nFloats, 'Floats.'

    
end subroutine

Function Float_Misfit(ObsState)
    double precision, dimension(Float_nObs,nj_DAstate, ni_DAstate), intent(in) :: ObsState
    double precision, dimension(Float_nObs,nj_DAstate, ni_DAstate) :: Float_Misfit
    
    integer indexi, indexj
    
    Float_Misfit=0.0d0
    if (.not.UseFloat) return
    
    indexj=0
    do indexi=1, nFloats
        if (Floats(indexi)%i==0) cycle
        if (DAMask(1, Floats(indexi)%j, Floats(indexi)%i)==0) cycle
        Float_Misfit(indexj+1:indexj+Floats(indexi)%nLevels, Floats(indexi)%j, Floats(indexi)%i) = &
            log(Floats(indexi)%tracer) - ObsState(indexj+1:indexj+Floats(indexi)%nLevels, Floats(indexi)%j, Floats(indexi)%i)        
        indexj=indexj+Floats(indexi)%nLevels
    end do
    
end function

function Float_H(State)
    
    double precision, dimension(nk_DAstate,nj_DAstate,ni_DAstate,ntra_DAstate), intent(in) :: State
    double precision, dimension(Float_nObs, nj_DAstate,ni_DAstate) :: Float_H
    
    integer indexi, indexj
    
    Float_H=0.0d0
    if (.not.UseFloat) return
    
    indexj=0
    do indexi=1, nFloats
        if (Floats(indexi)%i==0) cycle
        if (DAMask(1, Floats(indexi)%j, Floats(indexi)%i)==0) cycle
        Float_H(indexj+1:indexj+Floats(indexi)%nLevels, Floats(indexi)%j, Floats(indexi)%i) = &            
            log( sum(exp(State(ceiling(Floats(indexi)%k),Floats(indexi)%j,Floats(indexi)%i,FloatVariablesIndex)),2) * &
            (Floats(indexi)%k-floor(Floats(indexi)%k))  + &
            sum(exp(State(floor(Floats(indexi)%k),Floats(indexi)%j,Floats(indexi)%i,FloatVariablesIndex)),2) * &
            (1.0d0-(Floats(indexi)%k-floor(Floats(indexi)%k))) )
        indexj=indexj+Floats(indexi)%nLevels
    end do
    
end function

function Float_SqrtR1(HLi)
    double precision, dimension(Float_nObs, nj_DAstate, ni_DAstate), intent(in) :: HLi
    double precision, dimension(Float_n_SqrtR1, nj_DAstate, ni_DAstate) :: Float_SqrtR1
    
    integer indexi, indexj
    
    Float_SqrtR1=0.0d0
    if (.not.UseFloat) return
    
    indexj=0
    do indexi=1, nFloats
        if (Floats(indexi)%i==0) cycle
        if (DAMask(1, Floats(indexi)%j, Floats(indexi)%i)==0) cycle
        Float_SqrtR1(indexj+1:indexj+Floats(indexi)%nLevels, Floats(indexi)%j, Floats(indexi)%i) = &
            HLi(indexj+1:indexj+Floats(indexi)%nLevels,Floats(indexi)%j,Floats(indexi)%i) * Floats(indexi)%std1
        indexj=indexj+Floats(indexi)%nLevels
    end do
    
end function

function IsaFloatVar(name)
    implicit none
    
    character(LEN=*), intent(in) :: name
    logical :: IsaFloatVar
    
    if ((name.eq."P1l").or.(name.eq."P2l").or.(name.eq."P3l").or.(name.eq."P4l")) then
        IsaFloatVar=.true.
    else
        IsaFloatVar=.false.
    end if
        
end function

LOGICAL FUNCTION IsaDAFloat(datestring)
            
    CHARACTER(LEN=17), INTENT(IN) :: datestring
    ! LOCAL
    INTEGER I

    IsaDAFloat = .false.
    
    DO I=1, FloatTimes%N
        if (datestring.eq.FloatTimes%TimeStrings(I)) then
            IsaDAFloat = .true.
            return
        endif
    ENDDO
    
end function

end module
