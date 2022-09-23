!> Maps hard_open boundaries

module hard_open_mod
    use myalloc, only: lwp, find_index_var
    use bc_mod
    use bc_aux_mod
!    use myalloc
    implicit none

    private

    type, extends(bc) :: hard_open

        character(len=3) :: m_name ! ex: 'ope'
        integer :: m_n_vars
        character(len=20), allocatable, dimension(:) :: m_var_names
        character(len=20), allocatable, dimension(:) :: m_var_names_data
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
        double precision :: m_damping_factor

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
        procedure :: apply_dirichlet
        procedure :: apply_nudging
        procedure :: apply_phys
        procedure :: fix_diagnostic_vars ! overridden only by hard open class
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
        allocate(self%m_missing_var_names_idx(self%m_n_missing_vars))

        idx = 0
        do i = 1, n_tracers
            if (all(self%m_var_names_idx /= i)) then
                idx = idx + 1
                self%m_missing_var_names_idx(idx) = i
            endif
        enddo

    end subroutine set_missing_tracers



    !> Allocates and sets to 0 the 3d buffer matrix
    subroutine allocate_buffer(self)

        use modul_param, only: jpj, jpi, jpk

        implicit none

        ! TO DO: to be removed. Find a way to enable both testing and production code.
        ! integer, parameter :: jpk = 125
        ! integer, parameter :: jpj = 380
        ! integer, parameter :: jpi = 1085

        class(hard_open), intent(inout) :: self

        allocate(self%m_buffer(jpk, jpj, jpi))
        self%m_buffer = 0.0d0

    end subroutine allocate_buffer



    !> Defines a lookup matrix for the coordinates of the hard_open points
    subroutine set_hard_open_points(self)

        use modul_param, only: jpj, jpi, jpk

        implicit none

        ! TO DO: to be removed. Find a way to enable both testing and production code.
        ! integer, parameter :: jpk = 125
        ! integer, parameter :: jpj = 380
        ! integer, parameter :: jpi = 1085

        class(hard_open), intent(inout) :: self
        integer(4), allocatable, dimension(:, :) :: hard_open_points_aux ! TO DO: better use a 'stack' data structure
        integer :: i, j, k, counter

        ! TO DO: implement one more dimension to account for different geometries for different variables
        call readnc_slice_double(self%get_file_by_index(1), trim(self%m_var_names_data(1)), self%m_buffer)

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
        if (self%m_size > 0) then
            self%m_hard_open_points(:, :) = hard_open_points_aux(:, 1:self%m_size)
        endif

        deallocate(hard_open_points_aux)

    end subroutine set_hard_open_points



    !> Defines a lookup matrix for the coordinates of the neighbors
    subroutine set_neighbors(self)

        class(hard_open), intent(inout) :: self
        integer :: i

        ! TO DO: set condition for size = 0
        ! allcoate lookup matrix
        allocate(self%m_neighbors(3, self%m_size))

        if (self%m_size > 0) then
            
            select case(self%m_geometry)

                case(0) ! northern boundary
                    
                    do i = 1, self%m_size
                        self%m_neighbors(1, i) = self%m_hard_open_points(1, i)
                        self%m_neighbors(2, i) = self%m_hard_open_points(2, i) - 1
                        self%m_neighbors(3, i) = self%m_hard_open_points(3, i)
                    enddo
                
                case(1) ! eastern boundary
                    
                    do i = 1, self%m_size
                        self%m_neighbors(1, i) = self%m_hard_open_points(1, i) - 1
                        self%m_neighbors(2, i) = self%m_hard_open_points(2, i)
                        self%m_neighbors(3, i) = self%m_hard_open_points(3, i)
                    enddo
                
                case(2) ! southern boundary
                    
                    do i = 1, self%m_size
                        self%m_neighbors(1, i) = self%m_hard_open_points(1, i)
                        self%m_neighbors(2, i) = self%m_hard_open_points(2, i) + 1
                        self%m_neighbors(3, i) = self%m_hard_open_points(3, i)
                    enddo
                
                case(3) ! western boundary
                    
                    do i = 1, self%m_size
                        self%m_neighbors(1, i) = self%m_hard_open_points(1, i) + 1
                        self%m_neighbors(2, i) = self%m_hard_open_points(2, i)
                        self%m_neighbors(3, i) = self%m_hard_open_points(3, i)
                    enddo

            end select

        endif

    end subroutine set_neighbors



    ! 'bc_name' is used just to avoid system used symbol 'name'

    !> Target constructor

    !> Allocates and Initializes all the members that are added to the base class.
    subroutine init_members(self, bc_name, namelist_file, n_tracers)

        class(hard_open), intent(inout) :: self
        character(len=3), intent(in) :: bc_name
        character(len=7), intent(in) :: namelist_file
        integer, intent(in) :: n_tracers

        integer :: n_vars
        character(len=20), allocatable, dimension(:) :: vars
        integer(4) :: geometry
        double precision :: damping_coeff
        integer, parameter :: file_unit = 101 ! 100 for data files, 101 for boundary namelist files
        integer :: i
        namelist /vars_dimension/ n_vars
        namelist /core/ vars, geometry, damping_coeff

        self%m_name = bc_name

        ! read vars dimension parameters from namelist file
        open(unit=file_unit, file=namelist_file)
        rewind(file_unit)
        read(file_unit, vars_dimension)

        self%m_n_vars = n_vars

        ! allocate local arrays
        allocate(vars(self%m_n_vars))

        ! allocate class members
        allocate(self%m_var_names(self%m_n_vars))
        allocate(self%m_var_names_data(self%m_n_vars))
        allocate(self%m_var_names_idx(self%m_n_vars))

        ! read core parameters from namelist file
        rewind(file_unit)
        read(file_unit, core)

        do i = 1, self%m_n_vars
            self%m_var_names(i) = vars(i)
            self%m_var_names_data(i) = trim(self%m_var_names(i)) ! changed conventions, now just same as var names
            self%m_var_names_idx(i) = find_index_var(self%m_var_names(i))
        enddo

        self%m_geometry = geometry
        self%m_damping_factor = 1.0d0 / damping_coeff

        ! call delegated constructor - related procedures
        call self%set_missing_tracers(n_tracers)
        call self%allocate_buffer()
        call self%set_hard_open_points()
        call self%set_neighbors()

        allocate(self%m_values_dtatrc(self%m_size, 2, self%m_n_vars)) ! TO DO: which shape?
        self%m_values_dtatrc(:, :, :) = huge(self%m_values_dtatrc(1, 1, 1))
        allocate(self%m_values(self%m_size, self%m_n_vars))
        self%m_values(:, :) = huge(self%m_values(1, 1))

        ! deallocation
        deallocate(vars)

        ! close file
        close(unit=file_unit)

    end subroutine init_members



    ! TO DO: check if it is true that the constructor has to be always overloaded
    ! TO DO: final version of the constructor should receive everything from a single namelist
    ! 'bc_name' is used just to avoid system used symbol 'name'

    !> Default constructor

    !> Calls bc default constructor and target constructor.
    type(hard_open) function hard_open_default(bc_name, namelist_file, filenames_list, n_tracers)

        character(len=3), intent(in) :: bc_name
        character(len=7), intent(in) :: namelist_file
        character(len=22), intent(in) :: filenames_list
        integer, intent(in) :: n_tracers

        ! parent class constructor
        hard_open_default%bc = bc(filenames_list)

        call hard_open_default%init_members(bc_name, namelist_file, n_tracers)

        ! write(*, *) 'INFO: successfully called hard_open default constructor'

    end function hard_open_default



    ! TO DO: check if it is true that the constructor has to be always overloaded
    ! TO DO: final version of the constructor should receive everything from a single namelist
    ! 'bc_name' is used just to avoid system used symbol 'name'

    !> Periodic constructor

    !> Calls bc periodic constructor and target constructor.
    type(hard_open) function hard_open_year(bc_name, namelist_file, filenames_list, n_tracers, start_time_string, end_time_string)

        character(len=3), intent(in) :: bc_name
        character(len=7), intent(in) :: namelist_file
        character(len=22), intent(in) :: filenames_list
        integer, intent(in) :: n_tracers
        character(len=17), intent(in) :: start_time_string
        character(len=17), intent(in) :: end_time_string

        ! parent class constructor
        hard_open_year%bc = bc(filenames_list, start_time_string, end_time_string)

        call hard_open_year%init_members(bc_name, namelist_file, n_tracers)

        ! write(*, *) 'INFO: successfully called hard_open year constructor'

    end function hard_open_year



    !> Overridden from bc
    subroutine load(self, idx)

        class(hard_open), intent(inout) :: self
        integer, intent(in) :: idx
        integer :: i, j

        if (self%m_size > 0) then

            if (self%const_data()) then

                ! directly operate on self%m_values
                do i = 1, self%m_n_vars
                    call readnc_slice_double(self%get_file_by_index(idx), trim(self%m_var_names_data(i)), self%m_buffer)
                    do j = 1, self%m_size
                        self%m_values(j, i) = self%m_buffer( &
                            self%m_hard_open_points(3, j), &
                            self%m_hard_open_points(2, j), &
                            self%m_hard_open_points(1, j) &
                        )
                    enddo
                enddo

            else
                
                do i = 1, self%m_n_vars
                    call readnc_slice_double(self%get_file_by_index(idx), trim(self%m_var_names_data(i)), self%m_buffer)
                    do j = 1, self%m_size
                        self%m_values_dtatrc(j, 2, i) = self%m_buffer( &
                            self%m_hard_open_points(3, j), &
                            self%m_hard_open_points(2, j), &
                            self%m_hard_open_points(1, j) &
                        )
                    enddo
                enddo

            endif

        endif

    end subroutine load



    !> Overridden from bc
    subroutine swap(self)

        class(hard_open), intent(inout) :: self
        integer :: i, j

        if (self%m_size > 0) then

            ! if const_data(), do nothing
            if (.not.(self%const_data())) then
                
                do i = 1, self%m_n_vars
                    do j = 1, self%m_size
                        self%m_values_dtatrc(j, 1, i) = self%m_values_dtatrc(j, 2, i)
                    enddo
                enddo
            
            endif

        endif

    end subroutine swap



    !> Overridden from bc
    subroutine actualize(self, weight)

        class(hard_open), intent(inout) :: self
        double precision, intent(in) :: weight
        integer :: i, j

        if (self%m_size > 0) then

            ! if const_data(), do nothing
            if (.not.(self%const_data())) then
                
                do i = 1, self%m_n_vars
                    do j = 1, self%m_size
                        self%m_values(j, i) = &
                            (1.0 - weight) * self%m_values_dtatrc(j, 1, i) + weight * self%m_values_dtatrc(j, 2, i)
                    enddo
                enddo

            endif

        endif

    end subroutine actualize



    !> Overridden from bc

    !> Force all the open boundary points to the data values.
    !! They are not directly set to that value, but forced to it through a sort of hard nudging.
    !! This is just a more appropriate way to assign the values at the boundary cells,
    !! and a nudging on an arbitrary large set of cells can of course still be applied in the usual way through the decorator.
    !! A damping factor which should be of the same order of magnitude of a timestep is passed to the constructor.
    subroutine apply(self, e3t, n_tracers, trb, tra)

        use modul_param, only: jpk, jpj, jpi

        implicit none

        ! TO DO: to be removed. Find a way to enable both testing and production code.
        ! integer, parameter :: jpk = 125
        ! integer, parameter :: jpj = 380
        ! integer, parameter :: jpi = 1085

        class(hard_open), intent(inout) :: self
        double precision, dimension(jpk, jpj, jpi), intent(in) :: e3t
        integer, intent(in) :: n_tracers
        double precision, dimension(jpk, jpj, jpi, n_tracers), intent(in) :: trb
        double precision, dimension(jpk, jpj, jpi, n_tracers), intent(inout) :: tra
        integer :: i, j, idx_tracer, idx_i, idx_j, idx_k, idx_i_neigh, idx_j_neigh, idx_k_neigh
        double precision :: z_tracer

        if (self%m_size > 0) then

            ! First loop: provided variables.
            ! Values are taken from the data input and set through a strong nudging with a very small damping factor
            ! (which anyway is set from the boundary namelist)
            do i = 1, self%m_n_vars

                idx_tracer = self%m_var_names_idx(i)

                do j = 1, self%m_size

                    idx_i = self%m_hard_open_points(1, j)
                    idx_j = self%m_hard_open_points(2, j)
                    idx_k = self%m_hard_open_points(3, j)

                    z_tracer = self%m_damping_factor * ( &
                        self%m_values(j, i) - &
                        trb(idx_k, idx_j, idx_i, idx_tracer) &
                        )
                    tra(idx_k, idx_j, idx_i, idx_tracer) = tra(idx_k, idx_j, idx_i, idx_tracer) + z_tracer

                enddo

            enddo

            ! Second loop: missing variables.
            ! Values are taken from the neighbour cells and set through a strong nudging with a very small damping factor
            ! (which is the same used before)
            do i = 1, self%m_n_missing_vars

                idx_tracer = self%m_missing_var_names_idx(i)

                do j = 1, self%m_size

                    idx_i = self%m_hard_open_points(1, j)
                    idx_j = self%m_hard_open_points(2, j)
                    idx_k = self%m_hard_open_points(3, j)

                    idx_i_neigh = self%m_neighbors(1, j)
                    idx_j_neigh = self%m_neighbors(2, j)
                    idx_k_neigh = self%m_neighbors(3, j)

                    z_tracer = self%m_damping_factor * ( &
                        trb(idx_k_neigh, idx_j_neigh, idx_i_neigh, idx_tracer) - &
                        trb(idx_k, idx_j, idx_i, idx_tracer) &
                        )
                    tra(idx_k, idx_j, idx_i, idx_tracer) = tra(idx_k, idx_j, idx_i, idx_tracer) + z_tracer

                enddo

            enddo

        endif

    end subroutine apply

    subroutine apply_dirichlet(self)
        use modul_param, only: jpk, jpj, jpi
        use myalloc

        implicit none

        class(hard_open), intent(inout) :: self
        integer :: i, j, idx_tracer, idx_i, idx_j, idx_k

        if (self%m_size > 0) then
            ! First loop: provided variables on boundary points
            do i = 1, self%m_n_vars

                idx_tracer = self%m_var_names_idx(i)

                do j = 1, self%m_size

                    idx_i = self%m_hard_open_points(1, j)
                    idx_j = self%m_hard_open_points(2, j)
                    idx_k = self%m_hard_open_points(3, j)

                    tra(idx_k, idx_j, idx_i, idx_tracer) = self%m_values(j, i)

                enddo
             ! 2nd loop: provided variables on neighbors
                do j = 1, self%m_size

                    idx_i = self%m_neighbors(1, i)
                    idx_j = self%m_neighbors(2, j)
                    idx_k = self%m_neighbors(3, j)

                    tra(idx_k, idx_j, idx_i, idx_tracer) = self%m_values(j, i)

                enddo

            enddo
           endif
    end subroutine apply_dirichlet



    !> Overridden from bc
    subroutine apply_nudging(self, e3t, n_tracers, rst_tracers, trb, tra)

        use modul_param, only: jpk, jpj, jpi

        implicit none

        ! TO DO: to be removed. Find a way to enable both testing and production code.
        ! integer, parameter :: jpk = 125
        ! integer, parameter :: jpj = 380
        ! integer, parameter :: jpi = 1085

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
    subroutine apply_phys(self, lon, sponge_t, sponge_vel)

        use modul_param, only: jpk, jpj, jpi

        implicit none

        ! TO DO: to be removed. Find a way to enable both testing and production code.
        ! integer, parameter :: jpk = 125
        ! integer, parameter :: jpj = 380
        ! integer, parameter :: jpi = 1085

        class(hard_open), intent(inout) :: self
        double precision, dimension(jpj, jpi), intent(in) :: lon
        double precision, dimension(jpj, jpi), intent(out) :: sponge_t
        double precision, dimension(jpk, jpj, jpi), intent(out) :: sponge_vel

        if (lwp) write(*, *) 'WARN: sponge_t and sponge_val are left untouched by hard_open class'

    end subroutine apply_phys



    !> Overridden only by hard_open class

    !> It is supposed to assign values to the two diagnostic variables matrixes
    !> in the open boundary cells, where by default bfm is not executed.
    !> Values are set to be equal to the values of the neighbor cells.
    subroutine fix_diagnostic_vars(self)

        use modul_param, only: jpk, jpj, jpi
        use modul_param, only: jptra_dia, jptra_dia_2d
        use myalloc, only: tra_dia, tra_dia_2d

        implicit none

        class(hard_open), intent(inout) :: self

        integer :: i, jn, idx_i, idx_j, idx_k, idx_i_neigh, idx_j_neigh, idx_k_neigh

        do i = 1, self%m_size

            idx_i = self%m_hard_open_points(1, i)
            idx_j = self%m_hard_open_points(2, i)
            idx_k = self%m_hard_open_points(3, i)

            idx_i_neigh = self%m_neighbors(1, i)
            idx_j_neigh = self%m_neighbors(2, i)
            idx_k_neigh = self%m_neighbors(3, i)

            ! Set contributes on open boundary points equal to those on neighbor points for 3d matrix
            do jn = 1, jptra_dia
                tra_dia(idx_k, idx_j, idx_i,jn) = tra_dia(idx_k_neigh, idx_j_neigh, idx_i_neigh, jn)
            enddo

            ! Set contributes on open boundary points equal to those on neighbor points for 2d matrix
            do jn = 1, jptra_dia_2d
                tra_dia_2d(jn, idx_j, idx_i) = tra_dia_2d(jn, idx_j_neigh, idx_i_neigh)
            enddo

        enddo

    end subroutine fix_diagnostic_vars



    !> Destructor
    subroutine hard_open_destructor(self)

        class(hard_open), intent(inout) :: self

        if (allocated(self%m_var_names)) then
            deallocate(self%m_var_names)
            ! write(*, *) 'INFO: m_var_names deallocated'
        endif

        if (allocated(self%m_var_names_data)) then
            deallocate(self%m_var_names_data)
            ! write(*, *) 'INFO: m_var_names_data deallocated'
        endif

        if (allocated(self%m_var_names_idx)) then
            deallocate(self%m_var_names_idx)
            ! write(*, *) 'INFO: m_var_names_idx deallocated'
        endif

        if (allocated(self%m_missing_var_names_idx)) then
            deallocate(self%m_missing_var_names_idx)
            ! write(*, *) 'INFO: m_missing_var_names_idx deallocated'
        endif

        if (allocated(self%m_buffer)) then
            deallocate(self%m_buffer)
            ! write(*, *) 'INFO: m_buffer deallocated'
        endif

        if (allocated(self%m_hard_open_points)) then
            deallocate(self%m_hard_open_points)
            ! write(*, *) 'INFO: m_hard_open_points deallocated'
        endif

        if (allocated(self%m_neighbors)) then
            deallocate(self%m_neighbors)
            ! write(*, *) 'INFO: m_neighbors deallocated'
        endif

        if (allocated(self%m_values_dtatrc)) then
            deallocate(self%m_values_dtatrc)
            ! write(*, *) 'INFO: m_values_dtatrc deallocated'
        endif

        if (allocated(self%m_values)) then
            deallocate(self%m_values)
            ! write(*, *) 'INFO: m_values deallocated'
        endif

        ! parent class destructor
        call self%bc_destructor()

    end subroutine hard_open_destructor



end module hard_open_mod
