! developed by Simone Spada (sspada@ogs.it) at OGS

module Ens_Routines
    
implicit none
    
contains
    
    subroutine Ens_init
        use mpi
        
        use myalloc, &
            only: lwp
        use ogstm_mpi_module, &
            only: glcomm
        use Ens_Mem, &
            only: EnsDebug, EnsShareRestart, EnsRank, EnsSize, &
                Ens_restart_prefix, Ens_restart_ens_prefix, &
                Ens_ave_freq_1_prefix, Ens_ave_freq_1_ens_prefix, Ens_ave_freq_2_prefix, Ens_ave_freq_2_ens_prefix, &
                Ens_allocate
        use Ens_IO, &
            only: Ens_Init_IO, Ens_Finalize_IO, Ens_trcrst
        use Ens_Custom, &
            only: Ens_Init_DA
        use Ens_Utilities, &
            only: int2str
            
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
        
#endif
            
        end if
    end subroutine
    
    subroutine Ens_Finalize
        use Ens_Mem, &
            only: Ens_deallocate
        use Ens_IO, &
            only:  Ens_Finalize_IO
        use Ens_Custom, &
            only: Ens_Finalize_DA
        
        call Ens_deallocate
        call Ens_Finalize_IO
        
#ifdef ExecEnsDA
        call Ens_Finalize_DA
        
#endif
        
    end subroutine

end module
