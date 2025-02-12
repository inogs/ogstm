!> Inherits from bc implementing a decorator pattern

!> This pattern is a general pattern in Object-Oriented programming.
!! Since nudging is inheriting from bc, it is a bc.
!! Furthermore, associating its pointer to an already instantiated bc object of any kind
!! (both of base class or any of the derived classes)
!! it also has a bc, i.e. it can refer inside its methods directly to that object,
!! decorating it with additional features.
module nudging_mod
    use myalloc, only: lwp, find_index_var
    use bc_mod
    use rivers_mod
    use sponge_mod
    use hard_open_mod
    use bc_aux_mod

    implicit none

    private

    type, extends(bc) :: nudging

        class(bc), pointer :: m_bc_no_nudging => null()
        character(len=11) :: m_data_file ! 11 chars in order to handle names like 'bounmask.nc'
        integer :: m_n_nudging_vars
        character(len=20), allocatable, dimension(:) :: m_nudging_vars
        character(len=20), allocatable, dimension(:) :: m_nudging_vars_rst
        integer(4), allocatable, dimension(:) :: m_nudging_vars_idx ! tra_matrix_gib
        double precision, allocatable, dimension(:, :, :, :) :: m_rst ! resto
        double precision, allocatable, dimension(:) :: m_rst_corr ! restocorr
        double precision, allocatable, dimension(:, :, :, :) :: m_rst_tracers ! restotr

    contains

        ! delegated constructor
        procedure init_members
        ! base class methods
        procedure :: get_file_by_index
        procedure :: get_prev_idx
        procedure :: get_next_idx
        procedure :: set_current_interval
        procedure :: get_interpolation_factor
        procedure :: load
        procedure :: swap
        procedure :: actualize
        procedure :: apply
        procedure :: apply_nudging
        procedure :: apply_phys
        ! destructor
        procedure :: nudging_destructor

    end type nudging

    interface nudging
        module procedure nudging_default
    end interface nudging

    public :: nudging

contains



    !> Target constructor

    !> Allocates and Initializes all the members that are added to the base class.
    subroutine init_members(self, bc_no_nudging, namelist_file, n_tracers)

        use modul_param, only: jpk, jpj, jpi
        use netcdf
        use bc_aux_mod

        implicit none

        class(nudging), intent(inout) :: self
        class(bc), pointer, intent(in) :: bc_no_nudging
        character(len=7), intent(in) :: namelist_file
        integer, intent(in) :: n_tracers

        integer :: n_vars
        character(len=11) :: data_file ! 11 chars in order to handle names like 'bounmask.nc'
        character(len=20), allocatable, dimension(:) :: vars
        double precision, allocatable, dimension(:) :: rst_corr
        integer, parameter :: file_unit = 101 ! 100 for data files, 101 for boundary namelist files

        ! TO DO: to be removed. Find a way to enable both testing and production code.
        ! integer, parameter :: jpk = 125
        ! integer, parameter :: jpj = 380
        ! integer, parameter :: jpi = 1085

        integer :: i
        namelist /nudging_vars_dimension/ n_vars
        namelist /nudging_core/ data_file, vars, rst_corr

        ! pointer to bc_no_nudging association
        self%m_bc_no_nudging => bc_no_nudging

        ! read vars dimension parameters from namelist file
        open(unit=file_unit, file=namelist_file)
        rewind(file_unit)
        read(file_unit, nudging_vars_dimension)

        self%m_n_nudging_vars = n_vars

        ! allocate local arrays
        allocate(vars(self%m_n_nudging_vars))
        allocate(rst_corr(self%m_n_nudging_vars))

        ! allocate class members
        allocate(self%m_nudging_vars(self%m_n_nudging_vars))
        allocate(self%m_nudging_vars_rst(self%m_n_nudging_vars))
        allocate(self%m_nudging_vars_idx(self%m_n_nudging_vars))
        allocate(self%m_rst(jpk, jpj, jpi, self%m_n_nudging_vars))
        self%m_rst = huge(self%m_rst(1, 1, 1, 1))
        allocate(self%m_rst_corr(self%m_n_nudging_vars))
        allocate(self%m_rst_tracers(jpk, jpj, jpi, n_tracers))
        self%m_rst_tracers(:, :, :, :) = 0.

        ! read core parameters from namelist file
        rewind(file_unit)
        read(file_unit, nudging_core)

        self%m_data_file = data_file

        do i = 1, self%m_n_nudging_vars
            self%m_nudging_vars(i) = vars(i)
            self%m_nudging_vars_rst(i) = 're'//trim(self%m_nudging_vars(i))
            call readnc_slice_float(self%m_data_file, trim(self%m_nudging_vars_rst(i)), self%m_rst(:, :, :, i), 0)
            self%m_nudging_vars_idx(i) = find_index_var(self%m_nudging_vars(i))
            self%m_rst_corr(i) = rst_corr(i)
            self%m_rst_tracers(:, :, :, self%m_nudging_vars_idx(i)) = self%m_rst_corr(i) * self%m_rst(:, :, :, i)
        enddo

        ! deallocation
        deallocate(vars)
        deallocate(rst_corr)

        ! close file
        close(unit=file_unit)

    end subroutine init_members



    !> Default constructor

    !> Calls bc empty constructor and target constructor.
    !! No other constructors are needed so far.
    type(nudging) function nudging_default(bc_no_nudging, namelist_file, n_tracers)

        class(bc), pointer, intent(in) :: bc_no_nudging
        character(len=7), intent(in) :: namelist_file
        integer, intent(in) :: n_tracers

        ! parent class constructor
        nudging_default%bc = bc()

        call nudging_default%init_members(bc_no_nudging, namelist_file, n_tracers)

        ! write(*, *) 'INFO: successfully called nudging default constructor'

    end function nudging_default



    !> Overridden from bc.

    !> Redirects to the decorated object method.
    character(len=100) function get_file_by_index(self, idx)

        class(nudging), intent(in) :: self
        integer, intent(in) :: idx

        get_file_by_index = self%m_bc_no_nudging%get_file_by_index(idx)
        ! write(*, *) 'INFO: called get_file_by_index from nudging decorator'

    end function get_file_by_index



    !> Overridden from bc.

    !> Redirects to the decorated object method.
    integer function get_prev_idx(self)
        class(nudging), intent(in) :: self
        get_prev_idx = self%m_bc_no_nudging%get_prev_idx()
        ! write(*, *) 'INFO: called get_prev_idx from nudging decorator'
    end function get_prev_idx



    !> Overridden from bc.

    !> Redirects to the decorated object method.
    integer function get_next_idx(self)
        class(nudging), intent(in) :: self
        get_next_idx = self%m_bc_no_nudging%get_next_idx()
        ! write(*, *) 'INFO: called get_next_idx from nudging decorator'
    end function get_next_idx



    !> Overridden from bc.

    !> Redirects to the decorated object method.
    subroutine set_current_interval(self, current_time_string, new_data)

        class(nudging), intent(inout) :: self
        character(len=17), intent(in) :: current_time_string
        logical, optional, intent(out) :: new_data

        if (present(new_data)) then
            call self%m_bc_no_nudging%set_current_interval(current_time_string, new_data)
        else
            call self%m_bc_no_nudging%set_current_interval(current_time_string)
        endif
        ! write(*, *) 'INFO: called set_current_interval from nudging decorator'

    end subroutine set_current_interval



    !> Overridden from bc.

    !> Redirects to the decorated object method.
    double precision function get_interpolation_factor(self, current_time_string, new_data)

        class(nudging), intent(inout) :: self
        character(len=17), intent(in) :: current_time_string
        logical, optional, intent(out) :: new_data

        if (present(new_data)) then
            get_interpolation_factor = self%m_bc_no_nudging%get_interpolation_factor(current_time_string, new_data)
        else
            get_interpolation_factor = self%m_bc_no_nudging%get_interpolation_factor(current_time_string)
        endif
        ! write(*, *) 'INFO: called get_interpolation_factor from nudging decorator'

    end function get_interpolation_factor



    !> Overridden from bc.

    !> Redirects to the decorated object method.
    subroutine load(self, idx)

        class(nudging), intent(inout) :: self
        integer, intent(in) :: idx

        call self%m_bc_no_nudging%load(idx)
        ! write(*, *) 'INFO: called load from nudging decorator'

    end subroutine load



    !> Overridden from bc.

    !> Redirects to the decorated object method.
    subroutine swap(self)

        class(nudging), intent(inout) :: self

        call self%m_bc_no_nudging%swap()
        ! write(*, *) 'INFO: called swap from nudging decorator'

    end subroutine swap



    !> Overridden from bc.

    !> Redirects to the decorated object method.
    subroutine actualize(self, weight)

        class(nudging), intent(inout) :: self
        double precision, intent(in) :: weight

        call self%m_bc_no_nudging%actualize(weight)
        ! write(*, *) 'INFO: called actualize from nudging decorator'

    end subroutine actualize



    !> Overridden from bc.

    !> Redirects to the decorated object apply_nudging method.
    subroutine apply(self, e3t, n_tracers, trb, tra)

        use modul_param, only: jpk, jpj, jpi

        implicit none

        ! TO DO: to be removed. Find a way to enable both testing and production code.
        ! integer, parameter :: jpk = 125
        ! integer, parameter :: jpj = 380
        ! integer, parameter :: jpi = 1085

        class(nudging), intent(inout) :: self
        double precision, dimension(jpk, jpj, jpi), intent(in) :: e3t
        integer, intent(in) :: n_tracers
        double precision, dimension(jpk, jpj, jpi, n_tracers), intent(in) :: trb
        double precision, dimension(jpk, jpj, jpi, n_tracers), intent(inout) :: tra
        logical, save :: first = .true.

        if (first) then
           first=.false.
           !$acc enter data create(self%m_rst_tracers)
           !$acc update device(self%m_rst_tracers)
        endif

        call self%m_bc_no_nudging%apply_nudging(e3t, n_tracers, self%m_rst_tracers, trb, tra)
        ! write(*, *) 'INFO: called apply from nudging decorator'

    end subroutine apply



    !> Overridden from bc.

    !> Actually it does not do anything,
    !! since there is no need to apply an additional nudging to a nudging object.
    subroutine apply_nudging(self, e3t, n_tracers, rst_tracers, trb, tra)

        use modul_param, only: jpk, jpj, jpi

        implicit none

        ! TO DO: to be removed. Find a way to enable both testing and production code.
        ! integer, parameter :: jpk = 125
        ! integer, parameter :: jpj = 380
        ! integer, parameter :: jpi = 1085

        class(nudging), intent(inout) :: self
        double precision, dimension(jpk, jpj, jpi), intent(in) :: e3t
        integer, intent(in) :: n_tracers
        double precision, dimension(jpk, jpj, jpi, n_tracers), intent(in) :: rst_tracers
        double precision, dimension(jpk, jpj, jpi, n_tracers), intent(in) :: trb
        double precision, dimension(jpk, jpj, jpi, n_tracers), intent(inout) :: tra

        write(*, *) 'WARN: nudging class does not implement this method'
        write(*, *) 'WARN: attempt to apply nudging to nudging decorator'

    end subroutine apply_nudging



    !> Overridden from bc.

    !> Redirects to the decorated object method.
    subroutine apply_phys(self, lon, sponge_t, sponge_vel)

        use modul_param, only: jpk, jpj, jpi

        implicit none

        ! TO DO: to be removed. Find a way to enable both testing and production code.
        ! integer, parameter :: jpk = 125
        ! integer, parameter :: jpj = 380
        ! integer, parameter :: jpi = 1085

        class(nudging), intent(inout) :: self
        double precision, dimension(jpj, jpi), intent(in) :: lon
        double precision, dimension(jpj, jpi), intent(out) :: sponge_t
        double precision, dimension(jpk, jpj, jpi), intent(out) :: sponge_vel

        call self%m_bc_no_nudging%apply_phys(lon, sponge_t, sponge_vel)
        ! write(*, *) 'INFO: called apply_phys from nudging decorator'

    end subroutine apply_phys



    !> Destructor
    subroutine nudging_destructor(self)

        class(nudging), intent(inout) :: self

        class(bc), pointer :: bc_ptr ! auxiliary pointer to guarantee polymorphism

        bc_ptr => self%m_bc_no_nudging

        ! First call bc_no_nudging destructor according to its type.
        ! Only derived types (but not nudging (yet)) are supported
        select type (bc_ptr)

            class is (rivers)

                call bc_ptr%rivers_destructor()

            class is (sponge)

                call bc_ptr%sponge_destructor()

            class is (hard_open)

                call bc_ptr%hard_open_destructor()

        end select

        ! Then deallocate and nullyfy bc_no_nudging
        ! WARN: following line has been commented due tu a bug of the Intel compiler
        deallocate(self%m_bc_no_nudging)
        nullify(self%m_bc_no_nudging)

        ! write(*, *) 'INFO: m_bc_no_nudging deassociated'

        if (allocated(self%m_nudging_vars)) then
            deallocate(self%m_nudging_vars)
            ! write(*, *) 'INFO: m_nudging_vars deallocated'
        endif

        if (allocated(self%m_nudging_vars_rst)) then
            deallocate(self%m_nudging_vars_rst)
            ! write(*, *) 'INFO: m_nudging_vars_rst deallocated'
        endif

        if (allocated(self%m_nudging_vars_idx)) then
            deallocate(self%m_nudging_vars_idx)
            ! write(*, *) 'INFO: m_nudging_vars_idx deallocated'
        endif

        if (allocated(self%m_rst)) then
            deallocate(self%m_rst)
            ! write(*, *) 'INFO: m_rst deallocated'
        endif

        if (allocated(self%m_rst_corr)) then
            deallocate(self%m_rst_corr)
            ! write(*, *) 'INFO: m_rst_corr deallocated'
        endif

        if (allocated(self%m_rst_tracers)) then
            deallocate(self%m_rst_tracers)
            ! write(*, *) 'INFO: m_rst_tracers deallocated'
        endif

        ! parent class destructor
        call self%bc_destructor()

    end subroutine nudging_destructor



end module nudging_mod
