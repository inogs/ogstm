! developed by Simone Spada (sspada@ogs.it) at OGS

module Ens_Mem
    
implicit none
    integer, parameter :: EnsRankZero=0
    double precision, parameter :: Ens_Miss_val=1.0d20
    integer :: EnsDebug, EnsShareRestart
    integer :: EnsComm, EnsRank, EnsSize
    logical :: EnsSaveEachRestart, EnsSaveMeanRestart, EnsSaveEachAve, EnsSaveMeanAve, EnsAveDouble
    character(Len=100) :: Ens_restart_prefix, Ens_restart_ens_prefix, &
        Ens_ave_freq_1_prefix, Ens_ave_freq_1_ens_prefix, Ens_ave_freq_2_prefix, Ens_ave_freq_2_ens_prefix, &
        Ens_flux_prefix, Ens_flux_ens_prefix
    double precision, pointer, contiguous, dimension(:) :: DAstate 
    integer :: n_DAstate 
    integer :: win_DAstate
    double precision, pointer, contiguous, dimension(:,:) :: gl_DAstate
    
contains
    
    subroutine Ens_allocate
            
        !double precision, dimension(:), POINTER, contiguous :: member_pointer
        double precision, dimension(:,:), POINTER, contiguous :: global_pointer
        
        call Ens_shared_alloc(n_DAstate, DAstate, global_pointer, win_DAstate)
        !n_DAstate=jpk*jpj*jpi*jptra
        !DAstate(1:jpk,1:jpj,1:jpi,1:jptra)=>member_pointer
        gl_DAstate(1:n_DAstate, 0:EnsSize-1)=>global_pointer
        DAstate  = huge(DAstate(1)) 
        
    end subroutine
    
    subroutine Ens_deallocate
        use mpi 
        
        integer ierror
        
        CALL MPI_Win_free(win_DAstate, ierror)
        
    end subroutine
    
    Subroutine Ens_shared_alloc(member_size, member_pointer, global_pointer, window)
        USE mpi
        USE, INTRINSIC :: ISO_C_BINDING
        
        USE ogstm_mpi_module, &
            ONLY: glcomm
        
        INTEGER, intent(in) :: member_size
        double precision, dimension(:), POINTER, contiguous, intent(out) :: member_pointer
        double precision, dimension(:,:), POINTER, contiguous, intent(out) :: global_pointer
        INTEGER, intent(out) :: window
        
        INTEGER :: ierror
        INTEGER :: double_size
        INTEGER(KIND=MPI_ADDRESS_KIND) :: window_buffer_size
        TYPE(C_PTR) :: C_pointer
        
        CALL MPI_Type_size(MPI_real8, double_size, ierror)
        window_buffer_size = double_size * member_size
    
        CALL MPI_Win_allocate_shared(window_buffer_size, double_size, MPI_INFO_NULL, EnsComm, C_pointer, window, ierror)
        if (ierror/=0) then
            WRITE(*, *) 'Error with MPI_Win_allocate_shared. Aborting.'
            call MPI_abort(glcomm, 1, ierror)
        end if
 
        CALL C_F_POINTER(C_pointer, member_pointer, SHAPE = (/member_size/))
    
        call MPI_Win_shared_query(window, 0, window_buffer_size, double_size, C_pointer, ierror)
        
        CALL C_F_POINTER(C_pointer, global_pointer, SHAPE = (/member_size, EnsSize/))
    end subroutine
    
end module
