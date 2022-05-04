! developed by Simone Spada (sspada@ogs.it) at OGS

module Ens_Utilities
    implicit none
contains
    function str2int(str)
        implicit none

        character(len=*), intent(in) :: str
        integer :: str2int, ierr

        read(str,*,iostat=ierr) str2int
        if (ierr/=0) error stop "not a valid string in str2int"
    end function

    function int2str(number, ndigits)
        implicit none

        integer, intent(in) :: number, ndigits
        character(len=ndigits) :: int2str
        integer :: ierr
        character(len=10) :: str

        write(str,"(I0)") ndigits
        write(int2str,"(I0."//trim(str)//")",iostat=ierr) number
        if (ierr/=0) error stop "not a valid number in int2str"
    end function
    
    subroutine Ens_ReduceMean(window, member_size, gl_pointer)
        use mpi
        
        use myalloc, &
            only: myrank
        use Ens_Mem, &
            only: EnsRankZero, Ens_Miss_val, &
                EnsDebug, EnsRank, EnsSize
        
        integer, intent(in) :: window, member_size
        double precision, pointer, dimension(:,:), intent(inout) :: gl_pointer !dimension(member_size, 0:EnsSize-1)
         
        integer :: ierror
        integer :: istart, istop, indexi
        
        if (EnsDebug>1) write(*,*) '1st fence. EnsRank: ', EnsRank, ', myrank: ', myrank
        CALL MPI_Win_fence(0, window, ierror)
        
            if (EnsDebug>1) write(*,*) '1st fenced. EnsRank: ', EnsRank, ', myrank: ', myrank
            istart = 1 + (EnsRank*member_size)/EnsSize
            istop = ((EnsRank+1)*member_size)/EnsSize
            
            !gl_pointer(istart:istop, EnsRankZero) = sum(gl_pointer(istart:istop, :), dim=2)/EnsSize
            do indexi=istart, istop
                if (gl_pointer(indexi, EnsRankZero)<Ens_Miss_val) gl_pointer(indexi, EnsRankZero)=sum(gl_pointer(indexi, :))/EnsSize
            end do
            if (EnsDebug>1) write(*,*) 'reduction computed: ', EnsRank, ', myrank: ', myrank
            
        CALL MPI_Win_fence(0, window, ierror)
        if (EnsDebug>1) write(*,*) '2st fenced. EnsRank: ', EnsRank, ', myrank: ', myrank
        
    end subroutine
    
end module
