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
        procedure :: load
        procedure :: swap
        procedure :: actualize
        procedure :: apply
        procedure :: apply_phys
        ! destructor
        procedure :: nudging_destructor

    end type nudging

    interface nudging
        module procedure nudging_default
    end interface nudging

    public :: nudging

contains



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
        ! integer, parameter :: jpk = 70
        ! integer, parameter :: jpj = 65
        ! integer, parameter :: jpi = 182

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



    ! Default constructor invokes base class empty constructor;
    ! no other constructors are needed so far
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

    end function nudging_default



    subroutine load(self, idx)

        class(nudging), intent(inout) :: self
        integer, intent(in) :: idx

        call self%m_bc_no_nudging%load(idx)
        write(*, *) 'INFO: called load from nudging decorator'

    end subroutine load



    subroutine swap(self)

        class(nudging), intent(inout) :: self

        call self%m_bc_no_nudging%swap()
        write(*, *) 'INFO: called swap from nudging decorator'

    end subroutine swap



    subroutine actualize(self, weight)

        class(nudging), intent(inout) :: self
        double precision, intent(in) :: weight

        call self%m_bc_no_nudging%actualize(weight)
        write(*, *) 'INFO: called actualize from nudging decorator'

    end subroutine actualize



    subroutine apply(self, e3t, n_tracers, rst_tracers, trb, tra)

        use modul_param, only: jpk, jpj, jpi

        implicit none

        ! TO DO: to be removed. Find a way to enable both testing and production code.
        ! integer, parameter :: jpk = 70
        ! integer, parameter :: jpj = 65
        ! integer, parameter :: jpi = 182

        class(nudging), intent(inout) :: self
        double precision, dimension(jpk, jpj, jpi), intent(in) :: e3t
        integer, intent(in) :: n_tracers
        double precision, dimension(jpk, jpj, jpi, n_tracers), intent(in) :: rst_tracers
        double precision, dimension(jpk, jpj, jpi, n_tracers), intent(in) :: trb
        double precision, dimension(jpk, jpj, jpi, n_tracers), intent(inout) :: tra

        call self%m_bc_no_nudging%apply(e3t, n_tracers, self%m_rst_tracers, trb, tra)
        write(*, *) 'INFO: called apply from nudging decorator'

    end subroutine apply



    subroutine apply_phys(self, lat, sponge_t, sponge_vel)

        use modul_param, only: jpk, jpj, jpi

        implicit none

        ! TO DO: to be removed. Find a way to enable both testing and production code.
        ! integer, parameter :: jpk = 70
        ! integer, parameter :: jpj = 65
        ! integer, parameter :: jpi = 182

        class(nudging), intent(inout) :: self
        double precision, dimension(jpj, jpi), intent(in) :: lat
        double precision, dimension(jpj, jpi), intent(out) :: sponge_t
        double precision, dimension(jpk, jpj, jpi), intent(out) :: sponge_vel

        call self%m_bc_no_nudging%apply_phys(lat, sponge_t, sponge_vel)
        write(*, *) 'INFO: called apply_phys from nudging decorator'

    end subroutine apply_phys



    subroutine nudging_destructor(self)

        class(nudging), intent(inout) :: self

        ! just deassociate pointer (destructor will be invoked outside)
        nullify(self%m_bc_no_nudging)
        write(*, *) 'INFO: m_bc_no_nudging deassociated'

        if (allocated(self%m_nudging_vars)) then
            deallocate(self%m_nudging_vars)
            write(*, *) 'INFO: m_nudging_vars deallocated'
        endif

        if (allocated(self%m_nudging_vars_rst)) then
            deallocate(self%m_nudging_vars_rst)
            write(*, *) 'INFO: m_nudging_vars_rst deallocated'
        endif

        if (allocated(self%m_nudging_vars_idx)) then
            deallocate(self%m_nudging_vars_idx)
            write(*, *) 'INFO: m_nudging_vars_idx deallocated'
        endif

        if (allocated(self%m_rst)) then
            deallocate(self%m_rst)
            write(*, *) 'INFO: m_rst deallocated'
        endif

        if (allocated(self%m_rst_corr)) then
            deallocate(self%m_rst_corr)
            write(*, *) 'INFO: m_rst_corr deallocated'
        endif

        if (allocated(self%m_rst_tracers)) then
            deallocate(self%m_rst_tracers)
            write(*, *) 'INFO: m_rst_tracers deallocated'
        endif

        ! parent class destructor
        call self%bc_destructor()

    end subroutine nudging_destructor



end module nudging_mod
