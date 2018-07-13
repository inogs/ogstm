module bc_data_mod
    
    use calendar

    implicit none

    private

    ! type bc_data is supposed to contain values assumed by one single variable at a time;
    ! each boundary will thus host an array of bc_data,
    ! the size of which is will be given by the total number of variables.
    type bc_data
        ! file list related members
        integer :: m_n_files
        character(len=24), allocatable, dimension(:) :: m_files
        character(len=17), allocatable, dimension(:) :: m_time_strings
        ! time list related members
        integer :: m_n_times
        double precision, allocatable, dimension(:) :: m_times
    contains
        procedure :: get_file_by_index
        final :: destructor
    end type bc_data

    interface bc_data
        module procedure bc_data_default
        module procedure bc_data_year
        ! (optional) module procedure bc_data_full (time explicitly given)
    end interface bc_data

    public :: bc_data

contains

    ! TO DO: input file should be read with namelist syntax.
    ! TO DO: number of characters of 'files_namelist' is hard-coded (18)
    ! First entry in 'files_namelist' must be the number of files in the list.
    ! Default constructor infers times from file names:
    ! each file(5:21) must correspond to the fully qualified timestamp ('19700101-00:00:00').
    type(bc_data) function bc_data_default(files_namelist)

        character(len=18), intent(in) :: files_namelist
        integer, parameter :: file_unit = 100
        integer :: i ! counter

        ! open file
        open(unit=file_unit, file=files_namelist)

        ! get number of files and allocate memory accordingly
        read(file_unit, *) bc_data_default%m_n_files
        allocate(bc_data_default%m_files(bc_data_default%m_n_files))
        allocate(bc_data_default%m_time_strings(bc_data_default%m_n_files))
        bc_data_default%m_n_times = bc_data_default%m_n_files
        allocate(bc_data_default%m_times(bc_data_default%m_n_times))

        ! get file names and infer times
        read(file_unit, *) bc_data_default%m_files

        do i = 1, bc_data_default%m_n_files
            bc_data_default%m_time_strings(i) = bc_data_default%m_files(i)(5:21)
            bc_data_default%m_times(i) = datestring2sec(bc_data_default%m_time_strings(i))
        enddo

        ! close file
        close(unit=file_unit)

    end function bc_data_default

    ! TO DO: input files should be read with namelist syntax.
    ! TO DO: number of characters of 'files_namelist_year' is hard-coded (23)
    ! First entry in 'files_namelist' must be the number of files in the list.
    ! Year constructor infers times from file names and simulation duration:
    ! each file(9:21) must correspond to the fully qualified yearly timestamp ('0101-00:00:00').
    type(bc_data) function bc_data_year(files_namelist, start_time_string, end_time_string)

        character(len=23), intent(in) :: files_namelist
        character(len=17), intent(in) :: start_time_string
        character(len=17), intent(in) :: end_time_string
        integer, parameter :: file_unit = 100
        integer :: start_year, end_year, year, n_years
        character(len=4) :: year_string
        integer :: offset
        integer :: i ! counter

        ! open file
        open(unit=file_unit, file=files_namelist)

        ! get number of files and allocate memory accordingly
        read(file_unit, *) bc_data_year%m_n_files
        allocate(bc_data_year%m_files(bc_data_year%m_n_files))
        allocate(bc_data_year%m_time_strings(bc_data_year%m_n_files))
        
        ! get number of years and allocate memory accordingly
        read(start_time_string, '(i4)') start_year
        read(end_time_string, '(i4)') end_year
        n_years = end_year - start_year + 1
        bc_data_year%m_n_times = bc_data_year%m_n_files * n_years
        allocate(bc_data_year%m_times(bc_data_year%m_n_times))

        ! get file names and time strings
        read(file_unit, *) bc_data_year%m_files
        bc_data_year%m_time_strings(:) = bc_data_year%m_files(:)(5:21)

        ! infer times
        offset = 0
        do year = start_year, end_year
            write(year_string, '(i4)') year
            do i = 1, bc_data_year%m_n_files
                bc_data_year%m_time_strings(i)(1:4) = year_string
                bc_data_year%m_times(offset + i) = datestring2sec(bc_data_year%m_time_strings(i))
            enddo
            offset = offset + bc_data_year%m_n_files
        enddo

        ! reset 'bc_data_year%m_time_strings(:)(1:4)' to neutral string
        bc_data_year%m_time_strings(:)(1:4) = 'yyyy'

        ! close file
        close(unit=file_unit)

    end function bc_data_year

    ! This is supposed to match the given time to the right file, also with yearly data
    character(len=24) function get_file_by_index(self, idx)

        class(bc_data), intent(in) :: self
        integer, intent(in) :: idx
        integer :: file_idx
        
        if (idx > self%m_n_times) then
            write(*, *) 'ERROR: index out of range'
            get_file_by_index(:) = ' '
        else
            file_idx = mod(idx - 1, self%m_n_files) + 1
            get_file_by_index = self%m_files(file_idx)
        endif

    end function get_file_by_index

    subroutine destructor(self)

        type(bc_data), intent(inout) :: self

        if (allocated(self%m_files)) then
            deallocate(self%m_files)
            write(*, *) 'INFO: m_files deallocated'
        endif

        if (allocated(self%m_time_strings)) then
            deallocate(self%m_time_strings)
            write(*, *) 'INFO: m_time_strings deallocated'
        endif

        if (allocated(self%m_times)) then
            deallocate(self%m_times)
            write(*, *) 'INFO: m_times deallocated'
        endif

    end subroutine destructor

end module bc_data_mod
