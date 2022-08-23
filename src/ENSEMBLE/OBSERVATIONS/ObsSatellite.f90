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
        only: EnsRankZero, &
            EnsRank, EnsSize, &
            UseLocalObsDumping, &
            satfile_suffix, satvarname, &
            Ens_shared_alloc, &
            SatMultError
    use Ens_Custom, &   
        only: nk_DAstate, nj_DAstate, ni_DAstate, ntra_DAstate, &
            DAMask, DAVariablesIndex
    use LocalOperations, &
        only: LocalSpace
    
    
    implicit none
    
    logical, parameter :: LocalObs=.true.
    integer, parameter :: Sat_nObs=1, Sat_n_SqrtR1=1
    double precision, parameter :: fillvalue999=-999.0d0, fillValue = 1.0d20
    
    double precision, dimension(:,:), allocatable :: Sat_Data, Sat_var_additive, Sat_var_log, Sat_std1_ji
    integer, dimension(4) :: SatVariablesIndex
    type(dump_container) :: satTimes
    
    type(LocalSpace) :: LocalPatch
    
    
contains
    Subroutine Sat_Init
    
        !double precision, dimension(:), POINTER, contiguous :: member_pointer
        double precision, dimension(:,:), POINTER, contiguous :: global_pointer
        
        integer :: indexi, indexj
        
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
        
        if (IsaDASat(DateString)) then
            MONTH=DateString(5:6)
            DAY  =DateString(1:8)
            
            SATFILE   = 'SATELLITE/' // DAY // trim(satfile_suffix)
            VARFILE   = 'DA_static_data/VAR_SAT/var2D.' // MONTH // '.nc'
            call readnc_slice_double_2d(trim(SATFILE),trim(satvarname), Sat_Data)
            call readnc_slice_double_2d(trim(VARFILE),'variance', Sat_var_additive)
            !Sat_std1_ji=1.0d0/sqrt(Sat_var_additive(nldj:nlej, nldi:nlei))
            !Sat_var_log=0.3d0**2
            Sat_var_log=log(1.0d0+SatMultError)**2
            
            if (LocalObs) then
            
                do indexi=1,ni_DAstate
                    do indexj=1,nj_DAstate
                        tempvalue=Sat_Data(nldj-1+indexj,nldi-1+indexi)
                        if ( (.not.(tempvalue.eq.tempvalue)).or.(tempvalue.eq.fillValue).or.(tempvalue.eq.fillvalue999).or.(DAMask(1, indexj, indexi)==0) ) then 
                            Sat_std1_ji(indexj,indexi)=0.0d0
                        else
                            Sat_std1_ji(indexj,indexi)=1.0d0
                        end if
                    end do
                end do
                
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
        end if
        
    end subroutine
    
    Function Sat_Misfit(ObsState)
        double precision, dimension(Sat_nObs,nj_DAstate, ni_DAstate), intent(in) :: ObsState
        double precision, dimension(Sat_nObs,nj_DAstate, ni_DAstate) :: Sat_Misfit
        
        integer indexi, indexj
        double precision tempvalue
        
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
