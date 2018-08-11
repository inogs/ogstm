module bc_mod

    use bc_data_mod

    implicit none

    private

    type bc
        ! Eventually this should became an array of m_bc_data, one for each variable
        type(bc_data), pointer :: m_bc_data => null()
    contains
        procedure :: get_file_by_index
        ! TO DO: add all the procedures that are common to all the derived classes
        procedure :: bc_destructor
    end type bc

    interface bc
        module procedure bc_default
        module procedure bc_year
    end interface bc

    public :: bc

contains

    type(bc) function bc_default(files_namelist)

        character(len=18), intent(in) :: files_namelist

        allocate(bc_default%m_bc_data)
        bc_default%m_bc_data = bc_data(files_namelist)

    end function bc_default

    type(bc) function bc_year(files_namelist, start_time_string, end_time_string)

        character(len=23), intent(in) :: files_namelist
        character(len=17), intent(in) :: start_time_string
        character(len=17), intent(in) :: end_time_string

        allocate(bc_year%m_bc_data)
        bc_year%m_bc_data = bc_data(files_namelist, start_time_string, end_time_string)

    end function bc_year

    character(len=24) function get_file_by_index(self, idx)

        class(bc), intent(in) :: self
        integer, intent(in) :: idx

        get_file_by_index = self%m_bc_data%get_file_by_index(idx)

    end function get_file_by_index

    subroutine bc_destructor(self)

        class(bc), intent(inout) :: self

        ! First call bc_data destructor
        call self%m_bc_data%bc_data_destructor()
        
        deallocate(self%m_bc_data)
        write(*, *) 'INFO: m_bc_data deallocated'
        nullify(self%m_bc_data)
        write(*, *) 'INFO: m_bc_data deassociated'

    end subroutine bc_destructor

end module bc_mod
