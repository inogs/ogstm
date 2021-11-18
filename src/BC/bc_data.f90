!> This is a wrapper over the data structures which contains the full dataset needed by the boundary.

!> The dataset consists of a series of files referring to some specific times
!! distributed along the whole simulation period.
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
        logical :: m_const ! true when a single file is provided, only for periodic files
        character(len=100), allocatable, dimension(:) :: m_files
        character(len=17), allocatable, dimension(:) :: m_time_strings
        ! time list related members
        integer :: m_n_times
        double precision, allocatable, dimension(:) :: m_times
        ! interpolation variables
        integer :: m_prev_idx
        integer :: m_next_idx
        logical :: m_new_interval
    contains
        procedure :: get_file_by_index
        procedure :: const
        procedure :: get_prev_idx
        procedure :: get_next_idx
        procedure :: set_current_interval
        procedure :: get_interpolation_factor
        procedure :: new_interval
        procedure :: bc_data_destructor
    end type bc_data

    interface bc_data
        module procedure bc_data_empty
        module procedure bc_data_default
        module procedure bc_data_year
    end interface bc_data

    public :: bc_data

contains



    !> bc_data empty constructor
    type(bc_data) function bc_data_empty()

        bc_data_empty%m_n_files = 0
        bc_data_empty%m_const = .false. ! always false for empty constructor
        bc_data_empty%m_n_times = 0
        bc_data_empty%m_prev_idx = 0
        bc_data_empty%m_next_idx = 0
        bc_data_empty%m_new_interval = .false.

        ! write(*, *) 'INFO: successfully called bc_data empty constructor'

    end function bc_data_empty



    ! TO DO: input file should be read with namelist syntax.
    ! TO DO: number of characters of 'filenames_list' is hard-coded (22)
    ! First entry in 'filenames_list' must be the number of files in the list.
    ! Default constructor infers times from file names:
    ! each file(5:21) must correspond to the fully qualified timestamp ('19700101-00:00:00').

    !> bc_data default constructor

    !> Infers the time each file refers to from a time string contained in the file name itself.
    !! In this case the files and times lists will have the same length.
    type(bc_data) function bc_data_default(filenames_list)

        character(len=22), intent(in) :: filenames_list
        integer, parameter :: file_unit = 100
        integer :: i ! counter
        integer :: pos

        bc_data_default%m_const = .false. ! always false for default constructor

        ! open file
        open(unit=file_unit, file=filenames_list)

        ! get number of files and allocate memory accordingly
        read(file_unit, *) bc_data_default%m_n_files
        allocate(bc_data_default%m_files(bc_data_default%m_n_files))
        allocate(bc_data_default%m_time_strings(bc_data_default%m_n_files))
        bc_data_default%m_n_times = bc_data_default%m_n_files
        allocate(bc_data_default%m_times(bc_data_default%m_n_times))

        ! get file names and infer times
        read(file_unit, *) bc_data_default%m_files
        do i = 1, bc_data_default%m_n_files
            pos = 4
            bc_data_default%m_time_strings(i) = bc_data_default%m_files(i)(pos+4:pos+20)
            bc_data_default%m_times(i) = datestring2sec(bc_data_default%m_time_strings(i))
        enddo

        ! initialize interpolation variables
        bc_data_default%m_prev_idx = 1
        bc_data_default%m_next_idx = bc_data_default%m_prev_idx + 1
        bc_data_default%m_new_interval = .true.

        ! close file
        close(unit=file_unit)

        ! write(*, *) 'INFO: successfully called bc_data default constructor'

    end function bc_data_default



    ! TO DO: input files should be read with namelist syntax.
    ! TO DO: number of characters of 'filenames_list' is hard-coded (27)
    ! First entry in 'filenames_list' must be the number of files in the list.
    ! Year constructor infers times from file names and simulation duration:
    ! each file(9:21) must correspond to the fully qualified yearly timestamp ('0101-00:00:00').

    !> bc_data periodic constructor

    !> It is used to handle yearly periodic boundary conditions (e.g. climatological).
    !! Here only the constant part of the time string (i.e. month, day, hour etc.) is inferred from the files,
    !! and the list of times is computed and replicated for every simulation year, form start to end year.
    !! Therefore, this constructor accepts two arguments more (simulation start and end times).
    type(bc_data) function bc_data_year(filenames_list, start_time_string, end_time_string)

        character(len=22), intent(in) :: filenames_list
        character(len=17), intent(in) :: start_time_string
        character(len=17), intent(in) :: end_time_string
        integer, parameter :: file_unit = 100
        integer :: start_year, end_year, year, n_years
        character(len=4) :: year_string
        integer :: offset
        integer :: i ! counter

        ! open file
        open(unit=file_unit, file=filenames_list)

        ! get number of files and allocate memory accordingly
        read(file_unit, *) bc_data_year%m_n_files
        allocate(bc_data_year%m_files(bc_data_year%m_n_files))
        allocate(bc_data_year%m_time_strings(bc_data_year%m_n_files))
        
        ! check wether periodic boundary data are constant or not
        if (bc_data_year%m_n_files == 1) then
            bc_data_year%m_const = .true.
        else
            bc_data_year%m_const = .false.
        endif

        ! get number of years and allocate memory accordingly
        read(start_time_string, '(i4)') start_year
        read(end_time_string, '(i4)') end_year
        start_year = start_year - 1 ! in order to provide also previous year files for interpolation, if necessary
        end_year = end_year + 1 ! in order to provide also next year files for interpolation, if necessary
        n_years = end_year - start_year + 1
        bc_data_year%m_n_times = bc_data_year%m_n_files * n_years
        allocate(bc_data_year%m_times(bc_data_year%m_n_times))

        ! get file names and time strings
        read(file_unit, *) bc_data_year%m_files
        bc_data_year%m_time_strings(:) = bc_data_year%m_files(:)(8:24)

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

        ! initialize interpolation variables
        bc_data_year%m_prev_idx = 1
        bc_data_year%m_next_idx = bc_data_year%m_prev_idx + 1
        bc_data_year%m_new_interval = .true.

        ! close file
        close(unit=file_unit)

        ! write(*, *) 'INFO: successfully called bc_data year constructor'

    end function bc_data_year



    !> Constant data getter

    !> Boolean getter to check whether data are constant or not
    logical function const(self)
        class(bc_data), intent(in) :: self
        const = self%m_const
    end function const



    ! This is supposed to match the given time to the right file, also with yearly data

    !> Getter for the data file given the time index
    character(len=100) function get_file_by_index(self, idx)

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



    !> Getter for the left extreme of the current time interval
    integer function get_prev_idx(self)
        class(bc_data), intent(in) :: self
        get_prev_idx = self%m_prev_idx
    end function get_prev_idx



    !> Getter for the right extreme of the current time interval
    integer function get_next_idx(self)
        class(bc_data), intent(in) :: self
        get_next_idx = self%m_next_idx
    end function get_next_idx



    !> Setter of the current time interval

    !> Just sets the current time interval (prev and next indexes) according to current_time_string,
    !! without computing the interpolation factor.
    subroutine set_current_interval(self, current_time_string)

        class(bc_data), intent(inout) :: self
        character(len=17), intent(in) :: current_time_string
        double precision :: current_time
        integer :: i

        if (self%m_const) then

            ! in this case, it does nothing, and acts like interval never changes
            ! prev and next idx are set respectively to 1 nad 2 by the constructor
            ! and are left untouched
            self%m_new_interval = .false.

        else
            
            current_time = datestring2sec(current_time_string)
            
            ! TO DO: provide some bound checks
            i = 1
            do while(self%m_times(i+1) < current_time)
                i = i + 1
            enddo

            if (i /= self%m_prev_idx) then
                self%m_prev_idx = i
                self%m_next_idx = self%m_prev_idx + 1
                self%m_new_interval = .true.
            else
                self%m_new_interval = .false.
            endif

        endif

    end subroutine set_current_interval



    !> Linear interpolation factor computation

    !> Computes and returns the linear interpolation factor,
    !! keeping track of the current time interval (prev and next indexes).
    double precision function get_interpolation_factor(self, current_time_string)

        class(bc_data), intent(inout) :: self
        character(len=17), intent(in) :: current_time_string
        double precision :: current_time, prev_time, next_time

        call self%set_current_interval(current_time_string)

        if (self%m_const) then

            ! no need either to set interval nor to compute factor for constant data
            get_interpolation_factor = 1.0d0

        else
    
            current_time = datestring2sec(current_time_string)
            prev_time = self%m_times(self%m_prev_idx)
            next_time = self%m_times(self%m_next_idx)
            
            get_interpolation_factor = (current_time - prev_time) / (next_time - prev_time)

        endif

    end function get_interpolation_factor



    !> New interval getter

    !> Boolean getter to check whether the interval has changed or not
    !! after the last call either to set_current_interval or to get_interpolation_factor.
    logical function new_interval(self)
        class(bc_data), intent(in) :: self
        new_interval = self%m_new_interval
    end function new_interval



    !> Destructor
    subroutine bc_data_destructor(self)

        class(bc_data), intent(inout) :: self

        if (allocated(self%m_files)) then
            deallocate(self%m_files)
            ! write(*, *) 'INFO: m_files deallocated'
        endif

        if (allocated(self%m_time_strings)) then
            deallocate(self%m_time_strings)
            ! write(*, *) 'INFO: m_time_strings deallocated'
        endif

        if (allocated(self%m_times)) then
            deallocate(self%m_times)
            ! write(*, *) 'INFO: m_times deallocated'
        endif

    end subroutine bc_data_destructor



end module bc_data_mod
