module yearly_boundary_file_module
    use boundary_file_module
    use calendar
    implicit none

    private

    type, extends(boundary_file) :: yearly_boundary_file
        char(len=13) :: m_part_timestamp ! year-independent timestamp, e.g.: '0101-00:00:00'
    contains
        procedure :: extend_time
    end type yearly_boundary_file

    interface yearly_boundary_file
        module procedure yearly_boundary_file_default
        module procedure yearly_boundary_file_w_time
    end interface yearly_boundary_file

    public :: yearly_boundary_file

contains

    ! default constructor infers year-independent timestamp from file name
    type(yearly_boundary_file) function yearly_boundary_file_default(file_name)
        char(len=24), intent(in) :: file_name
        yearly_boundary_file_default%m_file_name = file_name
        yearly_boundary_file_default%m_part_timestamp = yearly_boundary_file_default%m_file_name(9:21)
        yearly_boundary_file_default%m_time = -1.0 ! dummy value: time is known only after extend_time is called
    end function yearly_boundary_file_default

    type(yearly_boundary_file) function yearly_boundary_file_w_time(file_name, part_timestamp)
        char(len=24), intent(in) :: file_name
        char(len=13), intent(in) :: part_timestamp
        yearly_boundary_file_w_time%m_file_name = file_name
        yearly_boundary_file_w_time%m_part_timestamp = part_timestamp
        yearly_boundary_file_w_time%m_time = -1.0 ! dummy value: time is known only after extend_time is called
    end function yearly_boundary_file_default

    subroutine extend_time(self, timestamp)
        class(boundary_file), intent(inout) :: self
        char(len=17), intent(in) :: timestamp
        self%m_time = datestring2sec(timestamp(1:4)//self%m_part_timestamp)
    end subroutine extend_time
    
end module yearly_boundary_file_module
