module boundary_mod

    use bc_data_mod

    implicit none

    private

    type boundary
        type(boundary_data), pointer :: m_data
    contains
        ! procedure :: allocate_memory (domrea.f90:215-216,230:231): should be done in the constructor
        final :: destructor
    end type boundary

    interface boundary
        module procedure boundary_default
        module procedure boundary_yearly
    end interface boundary

    public :: boundary

contains

    type(boundary) function boundary_default(files_namelist)
        char(len=24), intent(in) :: files_namelist
        allocate(m_data)
        m_data = boundary_data(files_namelsit)
    end function boundary_default

    type(boundary) function boundary_year(files_namelist, start_time_string, end_time_string)
        char(len=24), intent(in) :: files_namelist
        char(len=17), intent(in) :: start_time_string
        char(len=17), intent(in) :: end_time_string
        allocate(m_data)
        m_data = boundary_data(files_namelsit, start_time_string, end_time_string)
    end function boundary_year

    ! subroutine allocate_memory(self)
    !     class(boundary), intent(in) :: self
    !     write(*, *) 'WARNING: base class (boundary) does not implement this method'
    ! end subroutine

    subroutine destructor(self)
        class(boundary), intent(in) :: self
        nullify(self%m_data)
        deallocate(self%m_data)
    end subroutine

end module boundary_mod
