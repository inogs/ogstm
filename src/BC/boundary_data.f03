module boundary_data_module
    
    use calendar

    implicit none

    private

    type boundary_data
        ! file list related members
        integer :: m_n_files
        char(len=24), allocatable, dimension(:) :: m_files
        char(len=17), allocatable, dimension(:) :: m_time_strings
        ! time list related members
        integer :: m_n_times
        double precision, allocatable, dimension(:) :: m_times
    contains
        final :: destructor
    end type boundary_data

    interface boundary_data
        module procedure boundary_data_default
        module procedure boundary_data_year
        ! (optional) module procedure boundary_data_full
    end interface boundary_data

    public :: boundary_data

contains

    ! TO DO: input file should be read with namelist syntax.
    ! 'files_namelist' must contain in the first row the number of files in the list.
    ! Default constructor infers times from file names:
    ! each file(5:21) must correspond to the fully qualified timestamp ('19700101-00:00:00').
    type(boundary_data) function boundary_data_default(files_namelsit)

        char(len=24), intent(in) :: files_namelist

        ! open file
        open(unit=0, file=files_namelist)

        ! get number of files and allocate memory accordingly
        read(0, *) boundary_data_default%m_n_files
        allocate(boundary_data_default%m_files(boundary_data_default%m_n_files))
        allocate(boundary_data_default%m_time_strings(boundary_data_default%m_n_files))
        boundary_data_default%m_n_times = boundary_data_default%m_n_files
        allocate(boundary_data_default%m_times(boundary_data_default%m_n_times))

        ! get file names and infer times
        read(0, *) boundary_data_default%m_files
        boundary_data_default%m_time_strings(:) = boundary_data_default%m_files(:)(5:21)
        boundary_data_default%m_times(:) = datestring2sec(boundary_data_default%m_time_strings(:))

        ! close file
        close(unit=0)

    end function boundary_data_default

    ! TO DO: input files should be read with namelist syntax.
    ! 'files_namelist' must contain in the first row the number of files in the list.
    ! Year constructor infers times from file names and simulation duration:
    ! each file(9:21) must correspond to the fully qualified 'yearly' timestamp ('0101-00:00:00').
    type(boundary_data) function boundary_data_year(files_namelsit, start_time_string, end_time_string)

        char(len=24), intent(in) :: files_namelist
        char(len=17), intent(in) :: start_time_string
        char(len=17), intent(in) :: end_time_string
        integer :: start_year, end_year, year, n_years
        char(len=4) :: year_string
        integer :: offset = 1

        ! open file
        open(unit=0, file=files_namelist)

        ! get number of files and allocate memory accordingly
        read(0, *) boundary_data_year%m_n_files
        allocate(boundary_data_year%m_files(boundary_data_year%m_n_files))
        allocate(boundary_data_year%m_time_strings(boundary_data_year%m_n_files))
        
        ! get number of years and allocate memory accordingly
        read(start_time_string, '(i4)') start_year
        read(end_time_string, '(i4)') end_year
        n_years = end_year - start_year
        boundary_data_year%m_n_times = boundary_data_year%m_n_files * n_years
        allocate(boundary_data_year%m_times(boundary_data_year%m_n_times))

        ! get file names and time strings
        read(0, *) boundary_data_year%m_files
        boundary_data_year%m_time_strings(:) = boundary_data_year%m_files(:)(5:21)

        ! infer times
        do year = start_year, end_year
            write(year_string, '(i4)') year
            boundary_data_year%m_time_strings(:)(1:4) = year_string
            boundary_data_year%m_times(offset : offset + boundary_data_year%m_n_files - 1) = &
                datestring2sec(boundary_data_year%m_time_strings(:))
            offset += boundary_data_year%m_n_files
        enddo

        ! close file
        close(unit=0)

    end function boundary_data_default

    subroutine destructor(self)
        class(boundary_data), intent(inout) :: self
        deallocate(self%m_files)
        deallocate(self%m_time_strings)
        deallocate(self%m_times)
    end subroutine destructor

end module boundary_data_module
