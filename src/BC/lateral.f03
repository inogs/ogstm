module lateral_module
    use boundary_class
    implicit none

    private

    type, extends(boundary) :: lateral
        integer :: field_name
        ! - put here also the list of variables in domrea.F90:161-167,
        !   together with their number (7 is hard conded in BC_mem.f90:94)
        ! - put also file 'bounmask.nc'. Actually, this file is needed only in this class (see later)
        ! - decide wether to put here 'BC/GIB_'//TC_GIB%TimeStrings(1)//'.nc' (or the proper list of files) or to handle from parent
    contains
        procedure :: method_name
        ! - put here the procedures to read file 'bounmask.nc' (domrea.F90:169-180)
        !   'bounmask.nc' anyway should not be needed, and all the parameters should be initialized inside the class
        final :: destructor
    end type lateral

    interface lateral
        module procedure lateral_default
    end interface lateral

    public :: lateral

contains

    type(lateral) function lateral_default
    end function lateral_default

    subroutine method_name(self)
        class(lateral), intent(in) :: self
    end subroutine
    
    subroutine destructor(self)
        class(lateral), intent(in) :: self
    end subroutine

end module lateral_module
