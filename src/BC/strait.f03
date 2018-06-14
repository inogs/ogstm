module strait_module
    use boundary_class
    implicit none

    private

    type, extends(boundary) :: strait
        integer :: field_name
        ! - put here also the list of variables in domrea.F90:161-167,
        !   together with their number (7 is hard conded in BC_mem.f90:94)
        ! - put also file 'bounmask.nc'. Actually, this file is needed only in this class (see later)
    contains
        procedure :: method_name
        ! - put here the procedures to read file 'bounmask.nc' (domrea.F90:169-180)
        !   'bounmask.nc' anyway should not be needed, and all the parameters should be initialized inside the class
        final :: destructor
    end type strait

    interface strait
        module procedure strait_default
    end interface strait

    public :: strait

contains

    type(strait) function strait_default
    end function strait_default

    subroutine method_name(self)
        class(strait), intent(in) :: self
    end subroutine
    
    subroutine destructor(self)
        class(strait), intent(in) :: self
    end subroutine

end module strait_module
