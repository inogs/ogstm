module nudging_mod

    use bc_mod

    implicit none

    private

    type, extends(bc) :: nudging

        class(bc), pointer :: m_bc_no_nudging => null()
        integer :: m_decoration

    contains

        ! base class methods
        procedure :: load
        procedure :: swap
        procedure :: actualize
        ! destructor
        procedure :: nudging_destructor

    end type nudging

    interface nudging
        module procedure nudging_default
    end interface nudging

    public :: nudging

contains



    type(nudging) function nudging_default(filenames_list, bc_no_nudging, decoration)

        character(len=22), intent(in) :: filenames_list
        class(bc), target, intent(in) :: bc_no_nudging
        integer, intent(in) :: decoration

        ! parent class constructor
        nudging_default%bc = bc(filenames_list)

        ! pointer to bc_no_nudging association
        nudging_default%m_bc_no_nudging => bc_no_nudging

        ! other class members
        nudging_default%m_decoration = decoration

    end function nudging_default



    subroutine load(self, idx)

        class(nudging), intent(inout) :: self
        integer, intent(in) :: idx

        call self%m_bc_no_nudging%load(idx)
        write(*, *) 'INFO: called load from decorated class with decoration = ', self%m_decoration

    end subroutine load



    subroutine swap(self)

        class(nudging), intent(inout) :: self

        call self%m_bc_no_nudging%swap()
        write(*, *) 'INFO: called swap from decorated class with decoration = ', self%m_decoration

    end subroutine swap



    subroutine actualize(self, weight)

        class(nudging), intent(inout) :: self
        double precision, intent(in) :: weight

        call self%m_bc_no_nudging%actualize(weight)
        write(*, *) 'INFO: called actualize from decorated class with decoration = ', self%m_decoration

    end subroutine actualize



    subroutine nudging_destructor(self)

        class(nudging), intent(inout) :: self

        ! just deassociate pointer (destructor will be invoked outside)
        nullify(self%m_bc_no_nudging)
        write(*, *) 'INFO: m_bc_no_nudging deassociated'

        ! parent class destructor
        call self%bc_destructor()

    end subroutine nudging_destructor



end module nudging_mod
