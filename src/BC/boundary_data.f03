module boundary_data_module
    
    use calendar
    implicit none

    private

    type boundary_data
        char(len=24) :: m_file_name
        char(len=17) :: m_timestamp
        double precision m_time
    contains
        procedure :: get_time
        procedure :: extend_year
    end type boundary_data

    interface boundary_data
        module procedure boundary_data_default
        module procedure boundary_data_full
        module procedure boundary_data_part
    end interface boundary_data

    public :: boundary_data

contains

    ! Default constructor with auto = .true. infers time from file name:
    ! in this case, file_name must be a valid filename, i.e. either:
    ! - with fully qualified timestamp ('19700101-00:00:00');
    ! - with fully qualified timestamp without year characters ('yyyy0101-00:00:00').
    type(boundary_data) function boundary_data_default(file_name, auto)
        char(len=24), intent(in) :: file_name
        logical, intent(in), optional :: auto
        boundary_data_default%m_file_name = file_name
        boundary_data_default%m_timestamp = 'yyyy0101-00:00:00' ! dummy value
        boundary_data_default%m_time = -1.0 ! dummy value
        if (present(auto)) then
            if (auto) then
                boundary_data_default%m_timestamp = boundary_data_default%m_file_name(5:21)
                if (boundary_data_default%m_timestamp(1:4) \= 'yyyy') then
                    boundary_data_default%m_time = datestring2sec(boundary_data_default%m_timestamp)
                endif
            endif
        endif
    end function boundary_data_default

    ! timestamp must be either:
    ! - a fully qualified timestamp ('19700101-00:00:00');
    ! - a fully qualified timestamp without year characters ('yyyy0101-00:00:00').
    type(boundary_data) function boundary_data_full(file_name, timestamp)
        char(len=24), intent(in) :: file_name
        char(len=17), intent(in) :: timestamp
        boundary_data_full%m_file_name = file_name
        boundary_data_full%m_timestamp = timestamp
        if (boundary_data_full%m_timestamp(1:4) = 'yyyy') then
            boundary_data_full%m_time = -1.0 ! dummy value
        else
            boundary_data_full%m_time = datestring2sec(boundary_data_full%m_timestamp)
    end function boundary_data_full

    double precision function get_time(self)
        class(boundary_data), intent(in) :: self
        get_time = self%m_time
    end function get_time

    subroutine extend_year(self, timestamp)
        class(boundary_data), intent(inout) :: self
        char(len=17), intent(in) :: timestamp
        self%m_timestamp(1:4) = timestamp(1:4)
        self%m_time = datestring2sec(self%m_timestamp)
    end subroutine extend_year

end module boundary_data_module
