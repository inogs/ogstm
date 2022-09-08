!> Maps sponge boundaries

!> It maps a boundary which is forced to be closed even if in the OGCM is not.
!! It needs to know how to modify the velocities at the boundary;
!! in particular it should be able to set them to zero at the boundary
!! and to adapt them to the OGCM values according to a given function.
module sponge_mod

    use bc_mod
    use bc_aux_mod

    implicit none

    private

    type, extends(bc) :: sponge

        ! TO DO: review names
        character(len=3) :: m_name ! ex: 'gib'
        integer :: m_n_vars ! BC_mem.f90:94
        character(len=20), allocatable, dimension(:) :: m_var_names ! domrea.f90:161-167
        character(len=20), allocatable, dimension(:) :: m_var_names_data ! bc_gib.f90:113
        integer(4), allocatable, dimension(:) :: m_var_names_idx ! tra_matrix_gib
        double precision, allocatable, dimension(:, :, :) :: m_buffer ! replaces m_aux, now it is a 3D matrix
        integer(4) :: m_size ! BC_mem.f90:21
        integer(4), allocatable, dimension(:, :) :: m_sponge_points ! a lookup matrix
        double precision, allocatable, dimension(:, :, :) :: m_values_dtatrc ! TO DO: find better name
        double precision, allocatable, dimension(:, :) :: m_values
        ! sponge-related members
        double precision :: m_alpha
        double precision :: m_reduction_value_t
        double precision :: m_length

    contains

        ! target constructor - related procedures
        procedure :: allocate_buffer
        procedure :: set_sponge_points
        ! target constructor
        procedure :: init_members ! (memory allocation also in domrea.f90:215-216)
        ! base class methods
        procedure :: load
        procedure :: swap
        procedure :: actualize
        procedure :: apply
        procedure :: apply_nudging
        procedure :: apply_phys
        ! destructor
        procedure :: sponge_destructor

    end type sponge

    interface sponge
        module procedure sponge_default
        module procedure sponge_year
    end interface sponge

    public :: sponge

contains



    !> Allocates and sets to 0 the 3d buffer matrix
    subroutine allocate_buffer(self)

        use modul_param, only: jpj, jpi, jpk

        implicit none

        ! TO DO: to be removed. Find a way to enable both testing and production code.
        ! integer, parameter :: jpk = 125
        ! integer, parameter :: jpj = 380
        ! integer, parameter :: jpi = 1085

        class(sponge), intent(inout) :: self

        allocate(self%m_buffer(jpk, jpj, jpi))
        self%m_buffer = 0.0d0

    end subroutine allocate_buffer



    !> Defines a lookup matrix for the coordinates of the river points
    subroutine set_sponge_points(self)

        use modul_param, only: jpj, jpi, jpk

        implicit none

        ! TO DO: to be removed. Find a way to enable both testing and production code.
        ! integer, parameter :: jpk = 125
        ! integer, parameter :: jpj = 380
        ! integer, parameter :: jpi = 1085

        class(sponge), intent(inout) :: self
        integer(4), allocatable, dimension(:, :) :: sponge_points_aux ! TO DO: better use a 'stack' data structure
        integer :: i, j, k, counter

        ! TO DO: implement one more dimension to account for different geometries for different variables
        ! TO DO: single or double precision?
        call read_double(self%get_file_by_index(1), trim(self%m_var_names_data(1)),0, self%m_buffer)

        ! This is set to the largest dimension.
        ! Once computed, the real size will be used instead to allocate the right amount of memory
        allocate(sponge_points_aux(3, jpj*jpi*jpk))

        counter = 0
        do i = 1, jpi
            do j = 1, jpj
                do k = 1, jpk
                    ! upper / lower bounds have been decreased by one order of magnitude to avoid floating point issues
                    if ((self%m_buffer(k, j, i) < 1.0d19) .and. (self%m_buffer(k, j, i) > -1.0d-1)) then
                        counter = counter + 1
                        sponge_points_aux(1, counter) = i
                        sponge_points_aux(2, counter) = j
                        sponge_points_aux(3, counter) = k
                    endif
                enddo
            enddo
        enddo
        self%m_size = counter

        ! TO DO: set condition for size = 0
        ! allcoate lookup matrix
        allocate(self%m_sponge_points(3, self%m_size))
        ! copy
        if (self%m_size > 0) then
            self%m_sponge_points(:, :) = sponge_points_aux(:, 1:self%m_size)
        endif

        deallocate(sponge_points_aux)

    end subroutine set_sponge_points



    ! subroutine init_members(self, bc_name, bounmask, n_vars, vars)
    ! 'bc_name' is used just to avoid system used symbol 'name'

    !> Target constructor

    !> Allocates and Initializes all the members that are added to the base class.
    subroutine init_members(self, bc_name, namelist_file)

        class(sponge), intent(inout) :: self
        character(len=3), intent(in) :: bc_name
        character(len=7), intent(in) :: namelist_file

        integer :: n_vars
        character(len=20), allocatable, dimension(:) :: vars
        integer(4), allocatable, dimension(:) :: var_names_idx
        double precision :: alpha
        double precision :: reduction_value_t
        double precision :: length
        integer, parameter :: file_unit = 101 ! 100 for data files, 101 for boundary namelist files
        integer :: i
        namelist /vars_dimension/ n_vars
        namelist /core/ vars, var_names_idx, alpha, reduction_value_t, length

        self%m_name = bc_name

        ! read vars dimension parameters from namelist file
        open(unit=file_unit, file=namelist_file)
        rewind(file_unit)
        read(file_unit, vars_dimension)

        self%m_n_vars = n_vars

        ! allocate local arrays
        allocate(vars(self%m_n_vars))
        allocate(var_names_idx(self%m_n_vars))

        ! allocate class members
        allocate(self%m_var_names(self%m_n_vars))
        allocate(self%m_var_names_data(self%m_n_vars))
        allocate(self%m_var_names_idx(self%m_n_vars))

        ! read core parameters from namelist file
        rewind(file_unit)
        read(file_unit, core)

        do i = 1, self%m_n_vars
            self%m_var_names(i) = vars(i)
            self%m_var_names_data(i) = self%m_name//'_'//trim(self%m_var_names(i))
            self%m_var_names_idx(i) = var_names_idx(i)
        enddo

        ! call delegated constructor - related procedures
        call self%allocate_buffer()
        call self%set_sponge_points()

        allocate(self%m_values_dtatrc(self%m_size, 2, self%m_n_vars)) ! domrea.f90:216; TO DO: which shape?
        self%m_values_dtatrc(:, :, :) = huge(self%m_values_dtatrc(1, 1, 1)) ! domrea.f90:216
        allocate(self%m_values(self%m_size, self%m_n_vars)) ! domrea.f90:216
        self%m_values(:, :) = huge(self%m_values(1, 1)) ! domrea.f90:216

        self%m_alpha = alpha
        self%m_reduction_value_t = reduction_value_t
        self%m_length = length

        ! deallocation
        deallocate(vars)
        deallocate(var_names_idx)

        ! close file
        close(unit=file_unit)

    end subroutine init_members



    ! TO DO: check if it is true that the constructor has to be always overloaded
    ! TO DO: final version of the constructor should receive everything from a single namelist
    ! 'bc_name' is used just to avoid system used symbol 'name'

    !> Default constructor

    !> Calls bc default constructor and target constructor.
    type(sponge) function sponge_default(bc_name, namelist_file, filenames_list)

        character(len=3), intent(in) :: bc_name
        character(len=7), intent(in) :: namelist_file
        character(len=22), intent(in) :: filenames_list

        ! parent class constructor
        sponge_default%bc = bc(filenames_list)

        call sponge_default%init_members(bc_name, namelist_file)

        ! write(*, *) 'INFO: successfully called sponge default constructor'

    end function sponge_default



    ! TO DO: check if it is true that the constructor has to be always overloaded
    ! TO DO: final version of the constructor should receive everything from a single namelist
    ! 'bc_name' is used just to avoid system used symbol 'name'

    !> Periodic constructor

    !> Calls bc periodic constructor and target constructor.
    type(sponge) function sponge_year(bc_name, namelist_file, filenames_list, start_time_string, end_time_string)

        character(len=3), intent(in) :: bc_name
        character(len=7), intent(in) :: namelist_file
        character(len=22), intent(in) :: filenames_list
        character(len=17), intent(in) :: start_time_string
        character(len=17), intent(in) :: end_time_string

        ! parent class constructor
        sponge_year%bc = bc(filenames_list, start_time_string, end_time_string)

        call sponge_year%init_members(bc_name, namelist_file)

        ! write(*, *) 'INFO: successfully called sponge year constructor'

    end function sponge_year



    !> Overridden from bc
    subroutine load(self, idx)

        class(sponge), intent(inout) :: self
        integer, intent(in) :: idx
        integer :: i, j

        !if (self%m_size > 0) then

            if (self%const_data()) then

                ! directly operate on self%m_values
                do i = 1, self%m_n_vars
                    call read_double(self%get_file_by_index(idx), trim(self%m_var_names_data(i)),0, self%m_buffer)
                    do j = 1, self%m_size
                        self%m_values(j, i) = self%m_buffer( &
                            self%m_sponge_points(3, j), &
                            self%m_sponge_points(2, j), &
                            self%m_sponge_points(1, j) &
                        )
                    enddo
                enddo

            else
                
                do i = 1, self%m_n_vars
                    call read_double(self%get_file_by_index(idx), trim(self%m_var_names_data(i)),0, self%m_buffer)
                    do j = 1, self%m_size
                        self%m_values_dtatrc(j, 2, i) = self%m_buffer( &
                            self%m_sponge_points(3, j), &
                            self%m_sponge_points(2, j), &
                            self%m_sponge_points(1, j) &
                        )
                    enddo
                enddo

            endif

        !endif

    end subroutine load



    !> Overridden from bc
    subroutine swap(self)

        class(sponge), intent(inout) :: self
        integer :: i, j

        if (self%m_size > 0) then

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

        class(sponge), intent(inout) :: self
        double precision, intent(in) :: weight
        integer :: i, j

        if (self%m_size > 0) then

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

    !> Actually it does not do anything, since the chosen policy is to use this class always in combination with a nudging.
    subroutine apply(self, e3t, n_tracers, trb, tra)

        use modul_param, only: jpk, jpj, jpi

        implicit none

        ! TO DO: to be removed. Find a way to enable both testing and production code.
        ! integer, parameter :: jpk = 125
        ! integer, parameter :: jpj = 380
        ! integer, parameter :: jpi = 1085

        class(sponge), intent(inout) :: self
        double precision, dimension(jpk, jpj, jpi), intent(in) :: e3t
        integer, intent(in) :: n_tracers
        double precision, dimension(jpk, jpj, jpi, n_tracers), intent(in) :: trb
        double precision, dimension(jpk, jpj, jpi, n_tracers), intent(inout) :: tra

        write(*, *) 'WARN: sponge class does not implement this method'
        write(*, *) 'WARN: attempt to apply boundary conditions on sponge boundary without nudging'

    end subroutine apply



    !> Overridden from bc
    subroutine apply_nudging(self, e3t, n_tracers, rst_tracers, trb, tra)

        use modul_param, only: jpk, jpj, jpi

        implicit none

        ! TO DO: to be removed. Find a way to enable both testing and production code.
        ! integer, parameter :: jpk = 125
        ! integer, parameter :: jpj = 380
        ! integer, parameter :: jpi = 1085

        class(sponge), intent(inout) :: self
        double precision, dimension(jpk, jpj, jpi), intent(in) :: e3t
        integer, intent(in) :: n_tracers
        double precision, dimension(jpk, jpj, jpi, n_tracers), intent(in) :: rst_tracers
        double precision, dimension(jpk, jpj, jpi, n_tracers), intent(in) :: trb
        double precision, dimension(jpk, jpj, jpi, n_tracers), intent(inout) :: tra
        integer :: i, j, idx_tracer, idx_i, idx_j, idx_k
        double precision :: z_tracer

        if (self%m_size > 0) then
            do i = 1, self%m_n_vars
                idx_tracer = self%m_var_names_idx(i)
                do j = 1, self%m_size
                    idx_i = self%m_sponge_points(1, j)
                    idx_j = self%m_sponge_points(2, j)
                    idx_k = self%m_sponge_points(3, j)
                    z_tracer = rst_tracers(idx_k, idx_j, idx_i, idx_tracer) * &
                        (self%m_values(j, i) - trb(idx_k, idx_j, idx_i, idx_tracer))
                    tra(idx_k, idx_j, idx_i, idx_tracer) = tra(idx_k, idx_j, idx_i, idx_tracer) + z_tracer
                enddo
            enddo
        endif

    end subroutine apply_nudging



    !> Overridden from bc

    !> Provides the modified values to adjust both the velocities and the scalar fields of the OGCM.
    !! TO DO: provide more flexibility both
    !! in the choice of the 'sponge' function and
    !! in the choice of the interested area.
    !! So far, only a rectangular area in the west boundary has been implemented.
    !> WARNING: sponge_t and sponge_vel will be overwritten every time the method is called.
    subroutine apply_phys(self, lon, sponge_t, sponge_vel)

        use modul_param, only: jpk, jpj, jpi

        implicit none

        ! TO DO: to be removed. Find a way to enable both testing and production code.
        ! integer, parameter :: jpk = 125
        ! integer, parameter :: jpj = 380
        ! integer, parameter :: jpi = 1085

        class(sponge), intent(inout) :: self
        double precision, dimension(jpj, jpi), intent(in) :: lon ! glamt
        double precision, dimension(jpj, jpi), intent(out) :: sponge_t ! spongeT
        double precision, dimension(jpk, jpj, jpi), intent(out) :: sponge_vel ! spongeVel

        integer :: i, j
        double precision :: reduction_value_vel ! reduction_value (only for sponge_vel)

        do i = 1, jpi
            do j = 1, jpj
                if (lon(j, i) < self%m_length) then
                    sponge_t(j, i) = self%m_reduction_value_t
                    reduction_value_vel = exp( -self%m_alpha * ((lon(j, i) - self%m_length)**2) )
                    sponge_vel(:, j, i) = reduction_value_vel
                endif
            enddo
        enddo

    end subroutine apply_phys



    !> Destructor
    subroutine sponge_destructor(self)

        class(sponge), intent(inout) :: self

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

        if (allocated(self%m_buffer)) then
            deallocate(self%m_buffer)
            ! write(*, *) 'INFO: m_buffer deallocated'
        endif

        if (allocated(self%m_sponge_points)) then
            deallocate(self%m_sponge_points)
            ! write(*, *) 'INFO: m_sponge_points deallocated'
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

    end subroutine sponge_destructor



end module sponge_mod
