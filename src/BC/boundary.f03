module boundary_module
    implicit none

    private

    type boundary
        integer :: field_name
    contains
        procedure :: method_name
        final :: destructor
    end type boundary

    interface boundary
        module procedure boundary_default
    end interface boundary

    public :: boundary

contains

    type(boundary) function boundary_default
    end function boundary_default

    subroutine method_name(self)
        class(boundary), intent(in) :: self
    end subroutine
    
    subroutine destructor(self)
        class(boundary), intent(in) :: self
    end subroutine

end module boundary_module
