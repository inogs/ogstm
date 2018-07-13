module rivers_mod
    use boundary_mod
    implicit none

    private

    type, extends(edge) :: rivers
        integer :: field_name
    contains
        procedure :: method_name
        ! - read indexes (call to readnc_int_1d(), domrea.f90:226)
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

end module rivers_mod
