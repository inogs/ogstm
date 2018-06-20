module boundary_file_module
    use calendar
    implicit none

    private

    type boundary_file
        char(len=24) :: m_file_name
        double precision m_time
    contains
        procedure :: get_time
    end type boundary_file

    interface boundary_file
        module procedure boundary_file_default
        module procedure boundary_file_w_time
    end interface boundary_file

    public :: boundary_file

contains

    ! default constructor infers time from file name
    type(boundary_file) function boundary_file_default(file_name)
        char(len=24), intent(in) :: file_name
        boundary_file_default%m_file_name = file_name
        boundary_file_default%m_time = datestring2sec(boundary_file_default%m_file_name(5:21))
    end function boundary_file_default

    type(boundary_file) function boundary_file_w_time(file_name, timestamp)
        char(len=24), intent(in) :: file_name
        char(len=17), intent(in) :: timestamp
        boundary_file_w_time%m_file_name = file_name
        boundary_file_w_time%m_time = datestring2sec(timestamp)
    end function boundary_file_default

    double precision function get_time(self)
        class(boundary_file), intent(in) :: self
        get_time = self%m_time
    end subroutine

end module boundary_file_module
