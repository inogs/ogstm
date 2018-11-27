!> Maps hard_open boundaries

module hard_open_mod

    use bc_mod
    use bc_aux_mod

    implicit none

    private

    type, extends(bc) :: hard_open

        character(len=3) :: m_name ! ex: 'ope'
        integer :: m_n_vars
        character(len=3), allocatable, dimension(:) :: m_var_names
        character(len=7), allocatable, dimension(:) :: m_var_names_data
        integer(4), allocatable, dimension(:) :: m_var_names_idx
        integer :: m_n_missing_vars
        integer(4), allocatable, dimension(:) :: m_missing_var_names_idx
        double precision, allocatable, dimension(:, :, :) :: m_buffer ! replaces m_aux, now it is a 3D matrix
        integer(4) :: m_size
        integer(4), allocatable, dimension(:, :) :: m_hard_open_points ! a lookup matrix
        integer(4), allocatable, dimension(:, :) :: m_neighbors ! a lookup matrix for the neighbors
        double precision, allocatable, dimension(:, :, :) :: m_values_dtatrc
        double precision, allocatable, dimension(:, :) :: m_values
        integer(4) :: m_geometry ! 0=North; 1=East; 2=South; 3=West

    contains

        ! target constructor - related procedures
        procedure :: set_missing_tracers
        procedure :: allocate_buffer
        procedure :: set_hard_open_points
        procedure :: set_neighbors
        ! target constructor
        procedure :: init_members
        ! base class methods
        procedure :: load
        procedure :: swap
        procedure :: actualize
        procedure :: apply
        procedure :: apply_nudging
        procedure :: apply_phys
        ! new method for the hard open boundary
        procedure :: set_null_flux
        ! destructor
        procedure :: hard_open_destructor

    end type hard_open

    interface hard_open
        module procedure hard_open_default
        module procedure hard_open_year
    end interface hard_open

    public :: hard_open

contains



    !> Initializes the list of tracers which are not provided to the boundary

    !> Initializes the list of tracers (only their indexes) which are not provided to the boundary.
    !! The list of missing tracers is needed to know to which apply the no-flux condition.
    subroutine set_missing_tracers(self, n_tracers)

        class(hard_open), intent(inout) :: self
        integer, intent(in) :: n_tracers
        integer :: i, idx

        self%m_n_missing_vars = n_tracers - self%m_n_vars
        allocate(self%m_missing_vars(self%m_n_missing_vars))

        idx = 0
        do i = 1, n_tracers
            if (all(self%m_var_names_idx /= i)) then
                idx = idx + 1
                self%m_missing_vars(idx) = i
            endif
        enddo

    end subroutine set_missing_tracers



    !> Allocates and sets to 0 the 3d buffer matrix
    subroutine allocate_buffer(self)

        use modul_param, only: jpj, jpi, jpk

        implicit none

        ! TO DO: to be removed. Find a way to enable both testing and production code.
        ! integer, parameter :: jpj = 65
        ! integer, parameter :: jpi = 182
        ! integer, parameter :: jpk = 70

        class(hard_open), intent(inout) :: self

        allocate(self%m_buffer(jpk, jpj, jpi))
        self%m_buffer = 0.0d0

    end subroutine allocate_buffer



    !> Defines a lookup matrix for the coordinates of the hard_open points
    subroutine set_hard_open_points(self)

        use modul_param, only: jpj, jpi, jpk

        implicit none

        ! TO DO: to be removed. Find a way to enable both testing and production code.
        ! integer, parameter :: jpj = 65
        ! integer, parameter :: jpi = 182
        ! integer, parameter :: jpk = 70

        class(hard_open), intent(inout) :: self
        integer(4), allocatable, dimension(:, :) :: hard_open_points_aux ! TO DO: better use a 'stack' data structure
        integer :: i, j, k, counter

        ! TO DO: implement one more dimension to account for different geometries for different variables
        call readnc_slice_double(self%get_file_by_index(1), self%m_var_names_data(1), self%m_buffer)

        ! This is set to the largest dimension.
        ! Once computed, the real size will be used instead to allocate the right amount of memory
        allocate(hard_open_points_aux(3, jpj*jpi*jpk))

        counter = 0
        do i = 1, jpi
            do j = 1, jpj
                do k = 1, jpk
                    ! upper / lower bounds have been decreased by one order of magnitude to avoid floating point issues
                    ! sea points are set to 0, lower bound has been increased to 1.0d-6
                    if ((self%m_buffer(k, j, i) < 1.0d19) .and. (self%m_buffer(k, j, i) > 1.0d-6)) then
                        counter = counter + 1
                        hard_open_points_aux(1, counter) = i
                        hard_open_points_aux(2, counter) = j
                        hard_open_points_aux(3, counter) = k
                    endif
                enddo
            enddo
        enddo
        self%m_size = counter

        ! TO DO: set condition for size = 0
        ! allcoate lookup matrix
        allocate(self%m_hard_open_points(3, self%m_size))
        ! copy
        self%m_hard_open_points(:, :) = hard_open_points_aux(:, 1:self%m_size)

        deallocate(hard_open_points_aux)

    end subroutine set_hard_open_points



    !> Defines a lookup matrix for the coordinates of the neighbors
    subroutine set_neighbors(self)

        class(hard_open), intent(inout) :: self
        integer :: i

        ! TO DO: set condition for size = 0
        ! allcoate lookup matrix
        allocate(self%m_neighbors(3, self%m_size))

        select case(self%m_geometry)
        
            case(0) ! northern boundary

                do i = 1, self%m_size
                    self%m_neighbors(1, j) = self%m_hard_open_points(1, j)
                    self%m_neighbors(2, j) = self%m_hard_open_points(2, j) - 1
                    self%m_neighbors(3, j) = self%m_hard_open_points(3, j)
                enddo

            case(1) ! eastern boundary

                do i = 1, self%m_size
                    self%m_neighbors(1, j) = self%m_hard_open_points(1, j) - 1
                    self%m_neighbors(2, j) = self%m_hard_open_points(2, j)
                    self%m_neighbors(3, j) = self%m_hard_open_points(3, j)
                enddo

            case(2) ! southern boundary

                do i = 1, self%m_size
                    self%m_neighbors(1, j) = self%m_hard_open_points(1, j)
                    self%m_neighbors(2, j) = self%m_hard_open_points(2, j) + 1
                    self%m_neighbors(3, j) = self%m_hard_open_points(3, j)
                enddo

            case(3) ! western boundary

                do i = 1, self%m_size
                    self%m_neighbors(1, j) = self%m_hard_open_points(1, j) + 1
                    self%m_neighbors(2, j) = self%m_hard_open_points(2, j)
                    self%m_neighbors(3, j) = self%m_hard_open_points(3, j)
                enddo

        end select

    end subroutine set_neighbors



    ! 'bc_name' is used just to avoid system used symbol 'name'

    !> Target constructor

    !> Allocates and Initializes all the members that are added to the base class.
    subroutine init_members(self, bc_name, n_vars, vars, var_names_idx, n_tracers, geometry)

        class(hard_open), intent(inout) :: self
        character(len=3) :: bc_name
        integer, intent(in) :: n_vars
        character(len=23), intent(in) :: vars ! 'N1p N3n O2o N5s O3c O3h'; TO DO: more flexible
        integer(4), dimension(n_vars), intent(in) :: var_names_idx
        integer, intent(in) :: n_tracers
        integer(4) :: geometry
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
            self%m_var_names_data(i) = self%m_var_names(i) ! changed conventions, now just same as var names
            self%m_var_names_idx(i) = var_names_idx(i)
        enddo

        ! call delegated constructor - related procedures
        call self%set_missing_tracers(n_tracers)
        call self%allocate_buffer()
        call self%set_hard_open_points()
        call self%set_neighbors()

        allocate(self%m_values_dtatrc(self%m_size, 2, self%m_n_vars)) ! TO DO: which shape?
        self%m_values_dtatrc(:, :, :) = huge(self%m_values_dtatrc(1, 1, 1))
        allocate(self%m_values(self%m_size, self%m_n_vars))
        self%m_values(:, :) = huge(self%m_values(1, 1))

        self%m_geometry = geometry

    end subroutine init_members



    ! TO DO: check if it is true that the constructor has to be always overloaded
    ! TO DO: final version of the constructor should receive everything from a single namelist
    ! 'bc_name' is used just to avoid system used symbol 'name'

    !> Default constructor

    !> Calls bc default constructor and target constructor.
    type(hard_open) function hard_open_default(files_namelist, bc_name, n_vars, vars, var_names_idx, n_tracers, geometry)

        character(len=22), intent(in) :: files_namelist
        character(len=3) :: bc_name
        integer, intent(in) :: n_vars
        character(len=23), intent(in) :: vars
        integer(4), dimension(n_vars), intent(in) :: var_names_idx
        integer, intent(in) :: n_tracers
        integer(4) :: geometry

        ! parent class constructor
        hard_open_default%bc = bc(files_namelist)

        call hard_open_default%init_members(bc_name, n_vars, vars, var_names_idx, n_tracers, geometry)

        write(*, *) 'INFO: successfully called hard_open default constructor'

    end function hard_open_default



    ! TO DO: check if it is true that the constructor has to be always overloaded
    ! TO DO: final version of the constructor should receive everything from a single namelist
    ! 'bc_name' is used just to avoid system used symbol 'name'

    !> Periodic constructor

    !> Calls bc periodic constructor and target constructor.
    type(hard_open) function hard_open_year(files_namelist, bc_name, n_vars, vars, var_names_idx, n_tracers, geometry &
            start_time_string, end_time_string)

        character(len=27), intent(in) :: files_namelist
        character(len=3) :: bc_name
        integer, intent(in) :: n_vars
        character(len=27), intent(in) :: vars
        integer(4), dimension(n_vars), intent(in) :: var_names_idx
        integer, intent(in) :: n_tracers
        integer(4) :: geometry
        character(len=17), intent(in) :: start_time_string
        character(len=17), intent(in) :: end_time_string

        ! parent class constructor
        hard_open_year%bc = bc(files_namelist, start_time_string, end_time_string)

        call hard_open_year%init_members(bc_name, n_vars, vars, var_names_idx, n_tracers, geometry)

        write(*, *) 'INFO: successfully called hard_open year constructor'

    end function hard_open_year



    !> Overridden from bc
    subroutine load(self, idx)

        class(hard_open), intent(inout) :: self
        integer, intent(in) :: idx
        integer :: i, j

        do i = 1, self%m_n_vars
            call readnc_slice_double(self%get_file_by_index(idx), self%m_var_names_data(i), self%m_buffer)
            do j = 1, self%m_size
                self%m_values_dtatrc(j, 2, i) = self%m_buffer( &
                    self%m_hard_open_points(3, j), &
                    self%m_hard_open_points(2, j), &
                    self%m_hard_open_points(1, j) &
                )
            enddo
        enddo

    end subroutine load



    !> Overridden from bc
    subroutine swap(self)

        class(hard_open), intent(inout) :: self
        integer :: i, j

        do i = 1, self%m_n_vars
            do j = 1, self%m_size
                self%m_values_dtatrc(j, 1, i) = self%m_values_dtatrc(j, 2, i)
            enddo
        enddo

    end subroutine swap



    !> Overridden from bc
    subroutine actualize(self, weight)

        class(hard_open), intent(inout) :: self
        double precision, intent(in) :: weight
        integer :: i, j

        do i = 1, self%m_n_vars
            do j = 1, self%m_size
                self%m_values(j, i) = (1.0 - weight) * self%m_values_dtatrc(j, 1, i) + weight * self%m_values_dtatrc(j, 2, i)
            enddo
        enddo

    end subroutine actualize



    !> Overridden from bc

    !> Force all the open boundary points to the data values
    subroutine apply(self, e3t, n_tracers, trb, tra)

        use modul_param, only: jpk, jpj, jpi

        implicit none

        ! TO DO: to be removed. Find a way to enable both testing and production code.
        ! integer, parameter :: jpk = 70
        ! integer, parameter :: jpj = 65
        ! integer, parameter :: jpi = 182

        class(hard_open), intent(inout) :: self
        double precision, dimension(jpk, jpj, jpi), intent(in) :: e3t
        integer, intent(in) :: n_tracers
        double precision, dimension(jpk, jpj, jpi, n_tracers), intent(in) :: trb
        double precision, dimension(jpk, jpj, jpi, n_tracers), intent(inout) :: tra
        integer :: i, j, idx_tracer, idx_i, idx_j, idx_k

        if (self%m_size > 0) then
            do i = 1, self%m_n_vars
                idx_tracer = self%m_var_names_idx(i)
                do j = 1, self%m_size
                    tra( &
                        self%m_hard_open_points(3, j), &
                        self%m_hard_open_points(2, j), &
                        self%m_hard_open_points(1, j), &
                        idx_tracer &
                    ) = self%m_values(j, i)
                enddo
            enddo
        endif

    end subroutine apply



    !> Overridden from bc
    subroutine apply_nudging(self, e3t, n_tracers, rst_tracers, trb, tra)

        use modul_param, only: jpk, jpj, jpi

        implicit none

        ! TO DO: to be removed. Find a way to enable both testing and production code.
        ! integer, parameter :: jpk = 70
        ! integer, parameter :: jpj = 65
        ! integer, parameter :: jpi = 182

        class(hard_open), intent(inout) :: self
        double precision, dimension(jpk, jpj, jpi), intent(in) :: e3t
        integer, intent(in) :: n_tracers
        double precision, dimension(jpk, jpj, jpi, n_tracers), intent(in) :: rst_tracers
        double precision, dimension(jpk, jpj, jpi, n_tracers), intent(in) :: trb
        double precision, dimension(jpk, jpj, jpi, n_tracers), intent(inout) :: tra

        write(*, *) 'WARN: apply_nudging not implemented yet'

    end subroutine apply_nudging



    !> Overridden from bc

    !> Actually it does not do anything, since velocity fields from the OGCM are not modified
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

        write(*, *) 'WARN: sponge_t and sponge_val are left untouched by hard_open class'

    end subroutine apply_phys



    !> New method developed for the open boundary

    !> The aim is to set the concentrations of the missing tracers on the open boundary cells
    !! to the same values of the neighbor cells, in order to guarantee a null flux condition.
    !! This method is supposed to be called immediately before the call to the advection subroutine,
    !! in order to have no advection for the missing tracers at the boundary.
    subroutine set_null_flux(self, n_tracers, tra)

        use modul_param, only: jpk, jpj, jpi

        implicit none

        ! TO DO: to be removed. Find a way to enable both testing and production code.
        ! integer, parameter :: jpk = 70
        ! integer, parameter :: jpj = 65
        ! integer, parameter :: jpi = 182

        class(hard_open), intent(inout) :: self
        integer, intent(in) :: n_tracers
        double precision, dimension(jpk, jpj, jpi, n_tracers), intent(inout) :: tra
        integer :: i, j, idx_tracer

        if (self%m_size > 0) then
            do i = 1, self%m_n_missing_vars
                idx_tracer = self%m_missing_var_names_idx(i)
                do j = 1, self%m_size
                    tra( &
                        self%m_hard_open_points(3, j), &
                        self%m_hard_open_points(2, j), &
                        self%m_hard_open_points(1, j), &
                        idx_tracer &
                    ) = tra( &
                        self%m_neighbors(3, j), &
                        self%m_neighbors(2, j), &
                        self%m_neighbors(1, j), &
                        idx_tracer &
                    )
                enddo
            enddo
        endif

    end subroutine set_flux_null



    !> Destructor
    subroutine hard_open_destructor(self)

        class(hard_open), intent(inout) :: self

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

        if (allocated(self%m_missing_var_names_idx)) then
            deallocate(self%m_missing_var_names_idx)
            write(*, *) 'INFO: m_missing_var_names_idx deallocated'
        endif

        if (allocated(self%m_buffer)) then
            deallocate(self%m_buffer)
            write(*, *) 'INFO: m_buffer deallocated'
        endif

        if (allocated(self%m_hard_open_points)) then
            deallocate(self%m_hard_open_points)
            write(*, *) 'INFO: m_hard_open_points deallocated'
        endif

        if (allocated(self%m_neighbors)) then
            deallocate(self%m_neighbors)
            write(*, *) 'INFO: m_neighbors deallocated'
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

    end subroutine hard_open_destructor



end module hard_open_mod
