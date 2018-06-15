module rivers_module
    use boundary_class
    implicit none

    private

    type, extends(boundary) :: rivers
        integer :: field_name
        ! - decide wether to put here 'BC/TIN_'//TC_TIN%TimeStrings(1)//'.nc' (or the proper list of files) or to handle from parent
    contains
        procedure :: method_name
        final :: destructor
    end type rivers

    interface rivers
        module procedure rivers_default
    end interface rivers

    public :: rivers

contains

    type(rivers) function rivers_default
    end function rivers_default

    subroutine method_name(self)
        class(rivers), intent(in) :: self
    end subroutine
    
    subroutine destructor(self)
        class(rivers), intent(in) :: self
    end subroutine

end module rivers_module
