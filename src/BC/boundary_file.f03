module boundary_file_module
    implicit none

    private

    type boundary_file
        char(len=24) :: file_name
        double precision, dimension(:) :: timestamp
    contains
        procedure :: 
        final :: destructor
    end type boundary_file

    interface boundary_file
        module procedure boundary_file_default
    end interface boundary_file

    public :: boundary_file

contains

    type(boundary_file) function boundary_file_default
    end function boundary_file_default

    subroutine method_name(self)
        class(boundary_file), intent(in) :: self
    end subroutine
    
    subroutine destructor(self)
        class(boundary_file), intent(in) :: self
    end subroutine

end module boundary_file_module
