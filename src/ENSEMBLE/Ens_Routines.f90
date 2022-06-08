! developed by Simone Spada (sspada@ogs.it) at OGS

module Ens_Routines
    use mpi
    
    use myalloc, &
        only: lwp
    use ogstm_mpi_module, &
        only: glcomm
    use Ens_Mem, &
        only: EnsDebug, EnsShareRestart, EnsRank, EnsSize, &
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
        
    
implicit none
    
contains
    
    subroutine Ens_init
            
#ifdef ExecSEIK
        use ciccio, &
            only: createbase, readbase, createens
#endif
        
        integer ierr
        
        call Ens_allocate
        call Ens_Init_IO
        
        if (EnsSize<=1) then
        
            call Ens_trcrst(trim(Ens_restart_prefix), trim(Ens_ave_freq_1_prefix), trim(Ens_ave_freq_2_prefix))
            
        else
        
            Select case (EnsShareRestart)
        
                case (0) !each member has is own restart
                    call Ens_trcrst(trim(Ens_restart_ens_prefix)//int2str(EnsRank, 3), &
                        trim(Ens_ave_freq_1_ens_prefix)//int2str(EnsRank, 3), &
                        trim(Ens_ave_freq_2_ens_prefix)//int2str(EnsRank, 3))
#ifdef ExecSEIK
                    call createbase
#endif
                    
                case (1) !same initial restart from all ensemble members
                    call Ens_trcrst(trim(Ens_restart_prefix), trim(Ens_ave_freq_1_prefix), trim(Ens_ave_freq_2_prefix))
#ifdef ExecSEIK
                    call readbase
                    call createens
#endif

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
            
        end if
        
    end subroutine
    
    subroutine Ens_Finalize
        
        call Ens_deallocate
        call Ens_Finalize_IO
        
#ifdef ExecEnsDA
        if (EnsSize>1) then
            call Ens_Finalize_DA
            call Ens_Finalize_obs
            call Ens_Finalize_Filter
        end if
#endif
        
    end subroutine
    
    subroutine EnsDA(DateString)
        
        character(len=17), intent(in) :: DateString
        
        logical DoDA
        integer ierr
        
        DoDA=IsaDataAssimilation(DATEstring)
        if (IsaForecast(DATEstring).or.DoDA) then
            
            if (EnsDebug>1) then
                call mpi_barrier(glcomm, ierr)
                if (lwp) write(*,*) 'Starting EnsDA'
                DAparttime = MPI_WTIME()
            end if
            
            call Ens_state2DA
            
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
                call Ens_DA2state
                call Ens_SaveRestarts(Ens_forecast_prefix, Ens_forecast_ens_prefix, DATEstring)
            end if
            
            if (DoDA) then
                call Ens_LoadObs(DateString)
                
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
            
            call Ens_DA2state
            
            if (EnsDebug>1) then
                call mpi_barrier(glcomm, ierr)
                if (lwp) write(*,*) 'Ens_DA2state done'
            end if
            
            if (EnsSaveAfterAnalysis.and.DoDA) call Ens_SaveRestarts(Ens_analysis_prefix, Ens_analysis_ens_prefix, DATEstring)
            
            if (EnsDebug>0) then
                call mpi_barrier(glcomm, ierr)
                DAparttime = MPI_WTIME() - DAparttime
                DAtottime   = DAtottime + DAparttime
                if (lwp) write(*,*) 'End EnsDA. Time in seconds: ', DAparttime
            end if
            
        end if
    
    end subroutine

end module
