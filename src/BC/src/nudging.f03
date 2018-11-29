!> Inherits from bc implementing a decorator pattern

!> This pattern is a general pattern in Object-Oriented programming.
!! Since nudging is inheriting from bc, it is a bc.
!! Furthermore, associating its pointer to an already instantiated bc object of any kind
!! (both of base class or any of the derived classes)
!! it also has a bc, i.e. it can refer inside its methods directly to that object,
!! decorating it with additional features.
module nudging_mod

    use bc_mod

    implicit none

    private

    type, extends(bc) :: nudging

        class(bc), pointer :: m_bc_no_nudging => null()
        character(len=11) :: m_data_file ! 11 chars in order to handle names like 'bounmask.nc'
        integer :: m_n_nudging_vars
        character(len=3), allocatable, dimension(:) :: m_nudging_vars
        character(len=5), allocatable, dimension(:) :: m_nudging_vars_rst
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
    subroutine init_members(self, bc_no_nudging, data_file, n_vars, vars, vars_idx, rst_corr, n_tracers)

        use modul_param, only: jpk, jpj, jpi
        use netcdf
        use bc_aux_mod

        implicit none

        class(nudging), intent(inout) :: self
        class(bc), target, intent(in) :: bc_no_nudging
        character(len=11), intent(in) :: data_file ! 11 chars in order to handle names like 'bounmask.nc'
        integer, intent(in) :: n_vars
        character(len=27), intent(in) :: vars ! 'O2o N1p N3n N5s O3c O3h N6r'; TO DO: more flexible
        integer(4), dimension(n_vars), intent(in) :: vars_idx
        double precision, dimension(n_vars), intent(in) :: rst_corr
        integer, intent(in) :: n_tracers

        ! TO DO: to be removed. Find a way to enable both testing and production code.
        ! integer, parameter :: jpk = 125
        ! integer, parameter :: jpj = 380
        ! integer, parameter :: jpi = 1085

        integer :: i, start_idx, end_idx

        ! pointer to bc_no_nudging association
        self%m_bc_no_nudging => bc_no_nudging

        self%m_data_file = data_file
        self%m_n_nudging_vars = n_vars

        allocate(self%m_nudging_vars(self%m_n_nudging_vars))
        allocate(self%m_nudging_vars_rst(self%m_n_nudging_vars))
        allocate(self%m_nudging_vars_idx(self%m_n_nudging_vars))
        allocate(self%m_rst(jpk, jpj, jpi, self%m_n_nudging_vars))
        self%m_rst = huge(self%m_rst(1, 1, 1, 1))
        allocate(self%m_rst_corr(self%m_n_nudging_vars))
        allocate(self%m_rst_tracers(jpk, jpj, jpi, n_tracers))
        self%m_rst_tracers(:, :, :, :) = 0.
        do i = 1, self%m_n_nudging_vars
            end_idx = 4*i - 1
            start_idx = end_idx - 2
            self%m_nudging_vars(i) = vars(start_idx:end_idx)
            self%m_nudging_vars_rst(i) = 're'//self%m_nudging_vars(i)
            self%m_nudging_vars_idx(i) = vars_idx(i)
            call readnc_slice_float(self%m_data_file, self%m_nudging_vars_rst(i), self%m_rst(:, :, :, i), 0)
            self%m_rst_corr(i) = rst_corr(i)
            self%m_rst_tracers(:, :, :, self%m_nudging_vars_idx(i)) = self%m_rst_corr(i) * self%m_rst(:, :, :, i)
        enddo

    end subroutine init_members



    !> Default constructor

    !> Calls bc empty constructor and target constructor.
    !! No other constructors are needed so far.
    type(nudging) function nudging_default(bc_no_nudging, data_file, n_vars, vars, vars_idx, rst_corr, n_tracers)

        class(bc), target, intent(in) :: bc_no_nudging
        character(len=11), intent(in) :: data_file ! 11 chars in order to handle names like 'bounmask.nc'
        integer, intent(in) :: n_vars
        character(len=27), intent(in) :: vars ! 'O2o N1p N3n N5s O3c O3h N6r'; TO DO: more flexible
        integer(4), dimension(n_vars), intent(in) :: vars_idx
        double precision, dimension(n_vars), intent(in) :: rst_corr
        integer, intent(in) :: n_tracers

        ! parent class constructor
        nudging_default%bc = bc()

        call nudging_default%init_members(bc_no_nudging, data_file, n_vars, vars, vars_idx, rst_corr, n_tracers)

        ! write(*, *) 'INFO: successfully called nudging default constructor'

    end function nudging_default



    !> Overridden from bc.

    !> Redirects to the decorated object method.
    character(len=24) function get_file_by_index(self, idx)

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
    subroutine apply_phys(self, lat, sponge_t, sponge_vel)

        use modul_param, only: jpk, jpj, jpi

        implicit none

        ! TO DO: to be removed. Find a way to enable both testing and production code.
        ! integer, parameter :: jpk = 125
        ! integer, parameter :: jpj = 380
        ! integer, parameter :: jpi = 1085

        class(nudging), intent(inout) :: self
        double precision, dimension(jpj, jpi), intent(in) :: lat
        double precision, dimension(jpj, jpi), intent(out) :: sponge_t
        double precision, dimension(jpk, jpj, jpi), intent(out) :: sponge_vel

        call self%m_bc_no_nudging%apply_phys(lat, sponge_t, sponge_vel)
        ! write(*, *) 'INFO: called apply_phys from nudging decorator'

    end subroutine apply_phys



    !> Destructor
    subroutine nudging_destructor(self)

        class(nudging), intent(inout) :: self

        ! just deassociate pointer (destructor will be invoked outside)
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
