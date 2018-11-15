!> Maps rivers boundaries

!> It maps boundaries that in the boundary classification are defined as rivers.
!! This means that its data files are supposed to contain
!! the values of a potentially time dependent tracer flow at the boundaries,
!! i.e. discharges of an amount of tracers from the rivers to the ocean,
!! with a given seasonal variability.
module rivers_mod

    use bc_mod
    use bc_aux_mod

    implicit none

    private

    type, extends(bc) :: rivers

        ! TO DO: review names
        character(len=3) :: m_name ! ex: 'riv'
        integer :: m_n_vars ! BC_mem.f90:95
        character(len=3), allocatable, dimension(:) :: m_var_names
        character(len=7), allocatable, dimension(:) :: m_var_names_data ! bc_tin.f90:116
        integer(4), allocatable, dimension(:) :: m_var_names_idx ! tra_matrix_riv
        double precision, allocatable, dimension(:, :) :: m_buffer ! replaces m_aux, now it is a 2D matrix
        integer(4) :: m_size
        integer(4), allocatable, dimension(:,:) :: m_river_points ! a lookup matrix
        double precision, allocatable, dimension(:, :, :) :: m_values_dtatrc ! TO DO: find better name
        double precision, allocatable, dimension(:, :) :: m_values

    contains

        ! target constructor - related procedures
        procedure :: allocate_buffer
        procedure :: set_river_points
        ! target constructor
        procedure :: init_members
        ! base class methods
        procedure :: load
        procedure :: swap
        procedure :: actualize
        procedure :: apply
        procedure :: apply_nudging
        procedure :: apply_phys
        ! destructor
        procedure :: rivers_destructor

    end type rivers

    interface rivers
        module procedure rivers_default
        module procedure rivers_year
    end interface rivers

    public :: rivers

contains



    !> Allocates and initializes to 0 2d buffer matrix
    subroutine allocate_buffer(self)

        use modul_param, only: jpj, jpi

        implicit none

        ! TO DO: to be removed. Find a way to enable both testing and production code.
        ! integer, parameter :: jpj = 65
        ! integer, parameter :: jpi = 182

        class(rivers), intent(inout) :: self

        allocate(self%m_buffer(jpj, jpi))
        self%m_buffer = 0.0d0

    end subroutine allocate_buffer



    !> Defines a lookup matrix for the coordinates of the river points
    subroutine set_river_points(self)

        use modul_param, only: jpj, jpi

        implicit none

        ! TO DO: to be removed. Find a way to enable both testing and production code.
        ! integer, parameter :: jpj = 65
        ! integer, parameter :: jpi = 182

        class(rivers), intent(inout) :: self
        integer(4), allocatable, dimension(:,:) :: river_points_aux ! TO DO: better use a 'stack' data structure
        integer :: i, j, counter

        ! TO DO: implement one more dimension to account for different geometries for different variables
        call readnc_slice_float_2d(self%get_file_by_index(1), self%m_var_names_data(1), self%m_buffer, 0)

        ! This is set to the largest dimension.
        ! Once computed, the real size will be used instead to allocate the right amount of memory
        allocate(river_points_aux(2, jpj*jpi))

        counter = 0
        do i = 1, jpi
            do j = 1, jpj
                if (self%m_buffer(j, i) < 1.0d20) then
                    counter = counter + 1
                    river_points_aux(1, counter) = jpi
                    river_points_aux(2, counter) = jpj
                endif
            enddo
        enddo
        self%m_size = counter

        ! TO DO: set condition for size = 0
        ! allcoate aux matrix
        allocate(self%m_river_points(2, self%m_size))
        ! copy
        self%m_river_points(:, :) = river_points_aux(:, 1:self%m_size)

        deallocate(river_points_aux)

    end subroutine set_river_points



    ! 'bc_name' is used just to avoid system used symbol 'name'

    !> Target constructor

    !> Allocates and Initializes all the members that are added to the base class.
    subroutine init_members(self, bc_name, n_vars, vars, var_names_idx)

        class(rivers), intent(inout) :: self
        character(len=3) :: bc_name
        integer, intent(in) :: n_vars
        character(len=23), intent(in) :: vars ! 'N3n N1p N5s O3c O3h O2o'; TO DO: more flexible
        integer(4), dimension(n_vars), intent(in) :: var_names_idx
        integer :: i, start_idx, end_idx

        self%m_name = bc_name
        self%m_n_vars = n_vars

        allocate(self%m_var_names(self%m_n_vars))
        allocate(self%m_var_names_data(self%m_n_vars))
        allocate(self%m_var_names_idx(self%m_n_vars))

        do i = 1, self%m_n_vars
            end_idx = 4*i - 1
            start_idx = end_idx - 2
            self%m_var_names(i) = vars(start_idx:end_idx)
            self%m_var_names_data(i) = self%m_name//'_'//self%m_var_names(i)
            self%m_var_names_idx(i) = var_names_idx(i)
        enddo

        ! call delegated constructor - related procedures
        call self%allocate_buffer()
        call self%set_river_points()

        allocate(self%m_values_dtatrc(2, self%m_size, self%m_n_vars)) ! domrea.f90:216; TO DO: which shape?
        self%m_values_dtatrc(:, :, :) = huge(self%m_values_dtatrc(1, 1, 1)) ! domrea.f90:231
        allocate(self%m_values(self%m_size, self%m_n_vars)) ! domrea.f90:231
        self%m_values(:, :) = huge(self%m_values(1, 1)) ! domrea.f90:231

        write(*, *) 'INFO: successfully called init_members'

    end subroutine init_members



    ! TO DO: check if it is true that the constructor has to be always overloaded
    ! TO DO: final version of the constructor should receive everything from a single namelist
    ! 'bc_name' is used just to avoid system used symbol 'name'

    !> Default constructor

    !> Calls bc default constructor and target constructor.
    type(rivers) function rivers_default(files_namelist, bc_name, n_vars, vars, var_names_idx)

        character(len=22), intent(in) :: files_namelist
        character(len=3) :: bc_name
        integer, intent(in) :: n_vars
        character(len=23), intent(in) :: vars
        integer(4), dimension(n_vars), intent(in) :: var_names_idx

        ! parent class constructor
        rivers_default%bc = bc(files_namelist)

        call rivers_default%init_members(bc_name, n_vars, vars, var_names_idx)

        write(*, *) 'INFO: successfully called rivers default constructor'

    end function rivers_default



    ! TO DO: check if it is true that the constructor has to be always overloaded
    ! TO DO: final version of the constructor should receive everything from a single namelist
    ! 'bc_name' is used just to avoid system used symbol 'name'

    !> Periodic constructor

    !> Calls bc periodic constructor and target constructor.
    type(rivers) function rivers_year(files_namelist, bc_name, n_vars, vars, var_names_idx, start_time_string, end_time_string)

        character(len=27), intent(in) :: files_namelist
        character(len=3) :: bc_name
        integer, intent(in) :: n_vars
        character(len=23), intent(in) :: vars
        integer(4), dimension(n_vars), intent(in) :: var_names_idx
        character(len=17), intent(in) :: start_time_string
        character(len=17), intent(in) :: end_time_string

        ! parent class constructor
        rivers_year%bc = bc(files_namelist, start_time_string, end_time_string)

        call rivers_year%init_members(bc_name, n_vars, vars, var_names_idx)

        write(*, *) 'INFO: successfully called rivers year constructor'

    end function rivers_year



    !> Overridden from bc
    subroutine load(self, idx)

        class(rivers), intent(inout) :: self
        integer, intent(in) :: idx
        integer :: i, j

        do i = 1, self%m_n_vars
            call readnc_slice_float_2d(self%get_file_by_index(idx), self%m_var_names_data(i), self%m_buffer, 0)
            do j = 1, self%m_size
                self%m_values_dtatrc(2, j, i) = self%m_buffer(self%m_river_points(2, j), self%m_river_points(1, j))
            enddo
        enddo

    end subroutine load



    !> Overridden from bc
    subroutine swap(self)

        class(rivers), intent(inout) :: self
        integer :: i, j

        do i = 1, self%m_n_vars
            do j = 1, self%m_size
                self%m_values_dtatrc(1, j, i) = self%m_values_dtatrc(2, j, i)
            enddo
        enddo

    end subroutine swap



    !> Overridden from bc
    subroutine actualize(self, weight)

        class(rivers), intent(inout) :: self
        double precision, intent(in) :: weight
        integer :: i, j

        do i = 1, self%m_n_vars
            do j = 1, self%m_size
                self%m_values(j, i) = (1.0 - weight) * self%m_values_dtatrc(1, j, i) + weight * self%m_values_dtatrc(2, j, i)
            enddo
        enddo

    end subroutine actualize



    !> Overridden from bc
    subroutine apply(self, e3t, n_tracers, trb, tra)

        use modul_param, only: jpk, jpj, jpi

        implicit none

        ! TO DO: to be removed. Find a way to enable both testing and production code.
        ! integer, parameter :: jpk = 70
        ! integer, parameter :: jpj = 65
        ! integer, parameter :: jpi = 182        

        class(rivers), intent(inout) :: self
        double precision, dimension(jpk, jpj, jpi), intent(in) :: e3t
        integer, intent(in) :: n_tracers
        double precision, dimension(jpk, jpj, jpi, n_tracers), intent(in) :: trb
        double precision, dimension(jpk, jpj, jpi, n_tracers), intent(inout) :: tra
        integer :: i, j, idx_tracer, idx_i, idx_j

        if (self%m_size > 0) then
            do i = 1, self%m_n_vars
                idx_tracer = self%m_var_names_idx(i)
                do j = 1, self%m_size
                    idx_i = self%m_river_points(1, j)
                    idx_j = self%m_river_points(2, j)
                    tra(1, idx_j, idx_i, idx_tracer) = tra(1, idx_j, idx_i, idx_tracer) + &
                        self%m_values(j, i) / e3t(1, idx_j, idx_i)
                enddo
            enddo
        endif

    end subroutine apply



    !> Overridden from bc

    !> Actually it does not do anything, since so far the model does not need this feature for this type of boundary.
    subroutine apply_nudging(self, e3t, n_tracers, rst_tracers, trb, tra)

        use modul_param, only: jpk, jpj, jpi

        implicit none

        ! TO DO: to be removed. Find a way to enable both testing and production code.
        ! integer, parameter :: jpk = 70
        ! integer, parameter :: jpj = 65
        ! integer, parameter :: jpi = 182        

        class(rivers), intent(inout) :: self
        double precision, dimension(jpk, jpj, jpi), intent(in) :: e3t
        integer, intent(in) :: n_tracers
        double precision, dimension(jpk, jpj, jpi, n_tracers), intent(in) :: rst_tracers
        double precision, dimension(jpk, jpj, jpi, n_tracers), intent(in) :: trb
        double precision, dimension(jpk, jpj, jpi, n_tracers), intent(inout) :: tra
        integer :: i, j, idx_tracer, idx_i, idx_j

        write(*, *) 'WARN: rivers class does not implement this method'

    end subroutine apply_nudging



    !> Overridden from bc

    !> Actually it does not do anything,
    !! since rivers boundaries, by definition, are closed also in the OGCM;
    !! velocities are already set to zero and there is no point in modifying them at the boundaries.
    subroutine apply_phys(self, lat, sponge_t, sponge_vel)

        use modul_param, only: jpk, jpj, jpi

        implicit none

        ! TO DO: to be removed. Find a way to enable both testing and production code.
        ! integer, parameter :: jpk = 70
        ! integer, parameter :: jpj = 65
        ! integer, parameter :: jpi = 182

        class(rivers), intent(inout) :: self
        double precision, dimension(jpj, jpi), intent(in) :: lat
        double precision, dimension(jpj, jpi), intent(out) :: sponge_t
        double precision, dimension(jpk, jpj, jpi), intent(out) :: sponge_vel

        write(*, *) 'WARN: sponge_t and sponge_val are left untouched by rivers class'

    end subroutine apply_phys



    !> Destructor
    subroutine rivers_destructor(self)

        class(rivers), intent(inout) :: self

        if (allocated(self%m_var_names)) then
            deallocate(self%m_var_names)
            write(*, *) 'INFO: m_var_names deallocated'
        endif

        if (allocated(self%m_var_names_data)) then
            deallocate(self%m_var_names_data)
            write(*, *) 'INFO: m_var_names_data deallocated'
        endif

        if (allocated(self%m_var_names_idx)) then
            deallocate(self%m_var_names_idx)
            write(*, *) 'INFO: m_var_names_idx deallocated'
        endif

        if (allocated(self%m_buffer)) then
            deallocate(self%m_buffer)
            write(*, *) 'INFO: m_buffer deallocated'
        endif

        if (allocated(self%m_river_points)) then
            deallocate(self%m_river_points)
            write(*, *) 'INFO: m_river_points deallocated'
        endif

        if (allocated(self%m_values_dtatrc)) then
            deallocate(self%m_values_dtatrc)
            write(*, *) 'INFO: m_values_dtatrc deallocated'
        endif

        if (allocated(self%m_values)) then
            deallocate(self%m_values)
            write(*, *) 'INFO: m_values deallocated'
        endif

        ! parent class destructor
        call self%bc_destructor()

    end subroutine rivers_destructor



end module rivers_mod
