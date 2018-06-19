module boundary_file_module
    implicit none

    private

    type boundary_file
        char(len=24) :: file_name
        char(len=17) :: timestamp
        double precision time
    contains
        procedure :: get_time
        final :: destructor
    end type boundary_file

    interface boundary_file
        module procedure boundary_file_default
    end interface boundary_file

    public :: boundary_file

contains

    type(boundary_file) function boundary_file_default
    end function boundary_file_default

    double precision function get_time(self)
        class(boundary_file), intent(in) :: self
        get_time = self%time
    end subroutine
    
    subroutine destructor(self)
        class(boundary_file), intent(in) :: self
    end subroutine

end module boundary_file_module