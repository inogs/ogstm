!> Class bc_set represents a set of boundary conditions

!> It is supposed to be used to store and operate on all the boundary conditions imposed to the model.
!! It features an array of pointers to bc class objects, containing all the boundary conditions.
module bc_set_mod

    use bc_mod
    use rivers_mod
    use sponge_mod
    use hard_open_mod
    use nudging_mod
    use bc_handle_mod

    implicit none

    private

    type bc_container
        class(bc), pointer :: content => null()
    end type bc_container

    type bc_set
        integer :: m_n_bcs
        type(bc_container), allocatable, dimension(:) :: m_bcs
    contains
        !procedure :: update
        !procedure :: apply
        !procedure :: apply_phys
        procedure :: bc_set_destructor
    end type bc_set

    interface bc_set
        module procedure bc_set_default
    end interface bc_set

    public :: bc_set

contains



    type(bc_set) function bc_set_default(bcs_namelist)

        character(len=14), intent(in) :: bcs_namelist
        integer, parameter :: file_unit = 101 ! 100 for data files, 101 for boundary namelist files
        character(len=47) :: bc_string
        integer :: i

        ! open file
        open(unit=file_unit, file=bcs_namelist)

        ! get number of boundaries and allocate memory accordingly
        read(file_unit, *) bc_set_default%m_n_bcs
        allocate(bc_set_default%m_bcs(bc_set_default%m_n_bcs))

        ! get info for each boundary and initialize boundaries accordingly
        do i = 1, bc_set_default%m_n_bcs
            read(file_unit, *) bc_string
            write(*, *) 'INFO: initializing boundary with bc_string = ', bc_string
            bc_set_default%m_bcs(i)%content => bc_init(bc_string)
        enddo

        ! close file
        close(unit=file_unit)

    end function bc_set_default



    subroutine bc_set_destructor(self)

        class(bc_set), intent(inout) :: self
        integer :: i

        do i = 1, self%m_n_bcs
            call bc_destructor_wrapper(self%m_bcs(i)%content)
        enddo

        deallocate(self%m_bcs)

    end subroutine bc_set_destructor



end module bc_set_mod
