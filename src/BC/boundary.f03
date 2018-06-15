module boundary_module
    implicit none

    private

    type boundary
        integer :: field_name
        ! - decide wether to put the list of files here or in the child classes
    contains
        procedure :: method_name
        ! - read indexes (call to readnc_int_1d(), domrea.f90:204-207,226)
        ! - count in subdomain (domrea.f90:209-211,228)
        ! - allocate memory (domrea.f90:215-216,230:231)
        ! - re-indexing (domrea.f90:218,233, which should be moved away with the 3D index)
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
