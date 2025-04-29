! developed by Simone Spada (sspada@ogs.it) at OGS

module Ens_Routines
    use mpi
    
    use myalloc, &
        only: lwp
    use ogstm_mpi_module, &
        only: glcomm
    use Ens_Mem, &
        only: EnsRankZero, &
            EnsDebug, EnsRank, EnsSize, &
            UseParams, EnsShareRestart, EnsShareParams, &
            Ens_restart_prefix, Ens_restart_ens_prefix, &
            Ens_ave_freq_1_prefix, Ens_ave_freq_1_ens_prefix, Ens_ave_freq_2_prefix, Ens_ave_freq_2_ens_prefix, &
            EnsSaveAfterForecast, EnsSaveAfterAnalysis, &
            Ens_forecast_prefix, Ens_forecast_ens_prefix, Ens_analysis_prefix, Ens_analysis_ens_prefix, &
            IsaForecast, IsaDataAssimilation, &
            Ens_allocate, Ens_deallocate
    use Ens_IO, &
        only: Ens_Init_IO, Ens_Finalize_IO, Ens_trcrst, Ens_SaveRestarts
    use Ens_Custom, &
        only: Ens_state2DA, Ens_DA2state, &
            Ens_Init_DA, Ens_Finalize_DA
    use Ens_Utilities, &
        only: int2str
    use Ens_Obs, &
        only: nObs, &
            Ens_Init_Obs, Ens_Finalize_obs, Ens_LoadObs
    use Ens_Filter, &
        only: Ens_Init_Filter, Ens_Finalize_Filter, EnsForecast, EnsAnalysis
    use Timers, &
        only: DAparttime, DAtottime
    use Ens_Params, &
        only: Ens_Init_Params, Ens_Finalize_Params, Ens_ReadParams, Ens_WriteParams!, Ens_SetParams
    use Ens_ParallelWriter, &
        only: PW_Init, PW_Finalize
        
    
implicit none
    
contains
    
subroutine Ens_init
    
    integer ierr
    
    call Ens_allocate
    call Ens_Init_IO
    call PW_Init
    
#ifdef ExecEnsParams
    if (UseParams) call Ens_Init_Params
#endif

    if (EnsSize<=1) then
    
        call Ens_trcrst(trim(Ens_restart_prefix), trim(Ens_ave_freq_1_prefix), trim(Ens_ave_freq_2_prefix))
        
#ifdef ExecEnsParams
        if (UseParams) call Ens_ReadParams(trim(Ens_restart_prefix))
#endif
        
    else
    
        Select case (EnsShareRestart)
    
            case (0) !each member has its own restart
                call Ens_trcrst(trim(Ens_restart_ens_prefix)//int2str(EnsRank, 3), &
                    trim(Ens_ave_freq_1_ens_prefix)//int2str(EnsRank, 3), &
                    trim(Ens_ave_freq_2_ens_prefix)//int2str(EnsRank, 3))
                
            case (1) !same initial restart from all ensemble members
                call Ens_trcrst(trim(Ens_restart_prefix), trim(Ens_ave_freq_1_prefix), trim(Ens_ave_freq_2_prefix))

            case default
                if (lwp) write(*,*) "invalid EnsShareRestart value. Aborting."
                call mpi_barrier(glcomm,ierr)
                call MPI_abort(glcomm, 1, ierr)
                
        end select
        
#ifdef ExecEnsDA
        call Ens_Init_DA
        if (EnsDebug>1) then
            call mpi_barrier(glcomm,ierr)
            if (lwp) write(*,*) "Ens_Init_DA done!" 
        end if
        
        call Ens_Init_Obs
        if (EnsDebug>1) then
            call mpi_barrier(glcomm,ierr)
            if (lwp) write(*,*) "Ens_Init_Obs done!" 
        end if
        
        call Ens_Init_Filter
        if (EnsDebug>1) then
            call mpi_barrier(glcomm,ierr)
            if (lwp) write(*,*) "Ens_Init_Filter done!" 
        end if
#endif

#ifdef ExecEnsParams
        if (UseParams) then
            Select case (EnsShareParams)
        
                case (0) !each member has its own parameters
                    call Ens_ReadParams(trim(Ens_restart_ens_prefix)//int2str(EnsRank, 3), opt_subset_array=(/0/))
                    
                case (1) !same initial parameters from all ensemble members
                    call Ens_ReadParams(trim(Ens_restart_prefix), opt_subset_array=(/0/))

                case default
                    if (lwp) write(*,*) "invalid EnsShareParams value. Aborting."
                    call mpi_barrier(glcomm,ierr)
                    call MPI_abort(glcomm, 1, ierr)
                    
            end select
            
            call Ens_ReadParams(trim(Ens_restart_prefix), opt_subset_array=(/1/))
            
            if (EnsDebug>0) then
                call mpi_barrier(glcomm,ierr)
                if (lwp) write(*,*) "parameters loaded!" 
            end if
            
!                 call Ens_SetParams
        end if
#endif

        if (EnsDebug>0) then
            call mpi_barrier(glcomm,ierr)
            if (lwp) write(*,*) "Ens_Init done!" 
        end if
    end if
    
end subroutine

subroutine Ens_Finalize
    
    call Ens_deallocate
    call Ens_Finalize_IO
    call PW_Finalize
    
#ifdef ExecEnsDA
    if (EnsSize>1) then
        call Ens_Finalize_DA
        call Ens_Finalize_obs
        call Ens_Finalize_Filter
    end if
#endif

#ifdef ExecEnsParams
    if (UseParams) call Ens_Finalize_Params
#endif
    
end subroutine

subroutine EnsDA(DateString)
    
    character(len=17), intent(in) :: DateString
    
    logical DoDA
    integer ierr
    
!!!!debug lines
! if (UseParams) then
!     call Ens_WriteParams(DATEstring, trim(Ens_analysis_ens_prefix)//int2str(EnsRank, 3), opt_subset_array=(/0/))
!     
!     if (EnsRank==EnsRankZero) call Ens_WriteParams(DATEstring, trim(Ens_restart_prefix), opt_subset_array=(/1/))
!     
!     if (EnsDebug>0) then
!         call mpi_barrier(glcomm,ierr)
!         if (lwp) write(*,*) "parameters written!" 
!     end if
!     stop
! end if
!!!!!!!!!!!
    
    DoDA=IsaDataAssimilation(DATEstring)
    if (IsaForecast(DATEstring).or.DoDA) then
        
        if (EnsDebug>0) then
            call mpi_barrier(glcomm, ierr)
            if (lwp) write(*,*) 'Starting EnsDA'
            DAparttime = MPI_WTIME()
        end if
        
        call Ens_state2DA(DATEstring)
        
        if (EnsDebug>1) then
            call mpi_barrier(glcomm, ierr)
            if (lwp) write(*,*) 'Ens_state2DA done'
        end if
        
        call EnsForecast
        
        if (EnsDebug>1) then
            call mpi_barrier(glcomm, ierr)
            if (lwp) write(*,*) 'Forecast done'
        end if
        
        if (EnsSaveAfterForecast) then
            call Ens_DA2state(DATEstring)
            call Ens_SaveRestarts(Ens_forecast_prefix, Ens_forecast_ens_prefix, DATEstring)
        end if
        
        if (DoDA) then
            call Ens_LoadObs(DateString)
            
! if (EnsRank==EnsRankZero) call Ens_WriteParams(DATEstring, trim(Ens_analysis_prefix), opt_subset_array=(/1/))
! call mpi_barrier(glcomm,ierr)
! if (lwp) write(*,*) "parameters written!" 
! stop
            
            if (EnsDebug>1) then
                call mpi_barrier(glcomm, ierr)
                if (lwp) write(*,*) 'obs loaded'
            end if
            
            if (nObs==0) then
                if (lwp) write(*,*) DateString//': no observations available. Skipping DA.'
                DoDA=.false.
            else
                CALL EnsAnalysis
                
                if (EnsDebug>1) then
                    call mpi_barrier(glcomm, ierr)
                    if (lwp) write(*,*) 'Analysis done'
                end if
                
            end if
        endif
        
        call Ens_DA2state(DATEstring)
        
        if (EnsDebug>1) then
            call mpi_barrier(glcomm, ierr)
            if (lwp) write(*,*) 'Ens_DA2state done'
        end if
        
        if (EnsSaveAfterAnalysis.and.DoDA) then
        
            call Ens_SaveRestarts(Ens_analysis_prefix, Ens_analysis_ens_prefix, DATEstring)
        
            if (UseParams) then
                call Ens_WriteParams(DATEstring, trim(Ens_analysis_ens_prefix)//int2str(EnsRank, 3), opt_subset_array=(/0/))
                
                if (EnsRank==EnsRankZero) call Ens_WriteParams(DATEstring, trim(Ens_analysis_prefix), opt_subset_array=(/1/))
                
                if (EnsDebug>0) then
                    call mpi_barrier(glcomm,ierr)
                    if (lwp) write(*,*) "parameters written!" 
                end if
                
            end if
            
        end if
        
        if (EnsDebug>0) then
            call mpi_barrier(glcomm, ierr)
            DAparttime = MPI_WTIME() - DAparttime
            DAtottime   = DAtottime + DAparttime
            if (lwp) write(*,*) 'End EnsDA. Time in seconds: ', DAparttime
        end if
        
    end if

end subroutine

! subroutine test
!     
!     allocate(atest(3,EnsSize))
!     atest(:,1)=(/1,2,3/)
!     atest(:,2)=(/4,5,6/)
!     call test2(atest(1,:),1)
! 
! end subroutine
! 
! subroutine test2(a,n)
!     integer, intent(in) :: n
!     double precision, dimension(n,2), intent(in) :: a
!     integer ierr
!     write(*,*) a
!     call mpi_barrier(glcomm,ierr)
!     stop
! end subroutine

end module
