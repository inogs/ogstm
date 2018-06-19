module periodic_boundary_file_module
    implicit none

    private

    type periodic_boundary_file
        type(boundary_file) :: bf
    contains
        procedure :: get_time
        procedure :: extend_time
        final :: destructor
    end type periodic_boundary_file

    interface periodic_boundary_file
        module procedure periodic_boundary_file_default
    end interface periodic_boundary_file

    public :: periodic_boundary_file

contains

    type(periodic_boundary_file) function periodic_boundary_file_default
    end function periodic_boundary_file_default

    double precision function get_time(self)
        class(periodic_boundary_file), intent(in) :: self
        get_time = self%bf%get_time
    end subroutine

    subroutine extend_time(self, datestring)
    
    subroutine destructor(self)
        class(periodic_boundary_file), intent(in) :: self
    end subroutine

end module periodic_boundary_file_module
