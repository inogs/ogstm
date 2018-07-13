module bc_mod

    use bc_data_mod

    implicit none

    private

    type bc
        type(bc_data), pointer :: m_data
    contains
        ! procedure :: allocate_memory (domrea.f90:215-216,230:231): should be done in the constructor
        final :: destructor
    end type bc

    interface bc
        module procedure bc_default
        module procedure bc_yearly
    end interface bc

    public :: bc

contains

    type(bc) function bc_default(files_namelist)
        char(len=24), intent(in) :: files_namelist
        allocate(m_data)
        m_data = bc_data(files_namelsit)
    end function bc_default

    type(bc) function bc_year(files_namelist, start_time_string, end_time_string)
        char(len=24), intent(in) :: files_namelist
        char(len=17), intent(in) :: start_time_string
        char(len=17), intent(in) :: end_time_string
        allocate(m_data)
        m_data = bc_data(files_namelsit, start_time_string, end_time_string)
    end function bc_year

    ! subroutine allocate_memory(self)
    !     class(bc), intent(in) :: self
    !     write(*, *) 'WARNING: base class (bc) does not implement this method'
    ! end subroutine

    subroutine destructor(self)
        class(bc), intent(in) :: self
        nullify(self%m_data)
        deallocate(self%m_data)
    end subroutine

end module bc_mod
