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
        integer(4) :: m_global_size ! BC_mem.f90:20
        integer(4), allocatable, dimension(:) :: m_global_idxt ! BC_mem.f90:26; TO DO: find better name
        integer(4) :: m_size ! BC_mem.f90:21
        integer :: m_n_vars ! BC_mem.f90:94
        character(len=3), allocatable, dimension(:) :: m_var_names ! domrea.f90:161-167
        character(len=12), allocatable, dimension(:) :: m_var_names_idxt ! domrea.f90:204; TO DO: find better name
        character(len=7), allocatable, dimension(:) :: m_var_names_data ! bc_gib.f90:113
        integer(4), allocatable, dimension(:) :: m_var_names_idx ! tra_matrix_gib
        integer(4), allocatable, dimension(:, :) :: m_ridxt ! TO DO: find better name
        double precision, allocatable, dimension(:) :: m_aux
        double precision, allocatable, dimension(:, :, :) :: m_values_dtatrc ! TO DO: find better name
        double precision, allocatable, dimension(:, :) :: m_values
        ! sponge-related members
        double precision :: m_alpha
        double precision :: m_reduction_value_t
        double precision :: m_length

    contains

        ! delegated constructor - related procedures
        procedure :: set_global_size ! BC_mem.f90:70
        procedure :: set_global_idxt ! (call to readnc_int_1d(), domrea.f90:204-207)
        procedure :: set_size ! (domrea.f90:209-211)
        procedure :: reindex ! (domrea.f90:218, which should be moved away with the 3D index)
        ! delegated constructor
        procedure :: init_members ! (memory allocation also in domrea.f90:215-216)
        ! getters
        procedure :: get_global_size
        ! base class methods
        procedure :: load
        procedure :: swap
        procedure :: actualize
        procedure :: apply
        procedure :: apply_nudging
        procedure :: apply_phys ! TO DO: declare and implement also in base class
        ! destructor
        procedure :: sponge_destructor

    end type sponge

    interface sponge
        module procedure sponge_default
        module procedure sponge_year
    end interface sponge

    public :: sponge

contains



    !> Just a wrapper of 'getDimension'
    subroutine set_global_size(self)
        class(sponge), intent(inout) :: self
        call getDimension(self%get_file_by_index(1), self%m_var_names_idxt(1), self%m_global_size)
    end subroutine set_global_size



    !> Just a wrapper of 'readnc_int_1d'
    subroutine set_global_idxt(self)
        class(sponge), intent(inout) :: self
        allocate(self%m_global_idxt(self%m_global_size)) ! BC_mem.f90:111
        self%m_global_idxt(:) = huge(self%m_global_idxt(1))
        call readnc_int_1d(self%get_file_by_index(1), self%m_var_names_idxt(1), self%m_global_size, self%m_global_idxt)
    end subroutine set_global_idxt



    !> Just a wrapper of 'COUNT_InSubDomain_3d'
    subroutine set_size(self)
        class(sponge), intent(inout) :: self
        self%m_size = COUNT_InSubDomain_3d(self%m_global_size, self%m_global_idxt)
    end subroutine set_size



    !> Just a wrapper of 'RE_Indexing_3d'
    subroutine reindex(self)
        class(sponge), intent(inout) :: self
        call RE_Indexing_3d(self%m_global_size, self%m_global_idxt, self%m_size, self%m_ridxt)
    end subroutine reindex



    ! subroutine init_members(self, bc_name, bounmask, n_vars, vars)
    ! 'bc_name' is used just to avoid system used symbol 'name'

    !> Target constructor

    !> Allocates and Initializes all the members that are added to the base class.
    subroutine init_members(self, bc_name, n_vars, vars, var_names_idx, alpha, reduction_value_t, length)

        class(sponge), intent(inout) :: self
        character(len=3) :: bc_name
        integer, intent(in) :: n_vars
        character(len=27), intent(in) :: vars ! 'O2o N1p N3n N5s O3c O3h N6r'; TO DO: more flexible
        integer(4), dimension(n_vars), intent(in) :: var_names_idx
        double precision, intent(in) :: alpha
        double precision, intent(in) :: reduction_value_t
        double precision, intent(in) :: length
        integer :: i, start_idx, end_idx

        self%m_name = bc_name
        self%m_n_vars = n_vars

        allocate(self%m_var_names(self%m_n_vars))
        allocate(self%m_var_names_idxt(self%m_n_vars))
        allocate(self%m_var_names_data(self%m_n_vars))
        allocate(self%m_var_names_idx(self%m_n_vars))

        do i = 1, self%m_n_vars
            end_idx = 4*i - 1
            start_idx = end_idx - 2
            self%m_var_names(i) = vars(start_idx:end_idx)
            self%m_var_names_idxt(i) = self%m_name//'_idxt_'//self%m_var_names(i)
            self%m_var_names_data(i) = self%m_name//'_'//self%m_var_names(i)
            self%m_var_names_idx(i) = var_names_idx(i)
        enddo

        ! call delegated constructor - related procedures
        call self%set_global_size()
        call self%set_global_idxt()
        call self%set_size()

        allocate(self%m_ridxt(4, self%m_size)) ! domrea.f90:216
        self%m_ridxt(:, :) = huge(self%m_ridxt(1, 1)) ! domrea.f90:216
        call self%reindex() ! domrea.f90:218

        allocate(self%m_aux(self%m_global_size))
        self%m_aux(:) = huge(self%m_aux(1))
        allocate(self%m_values_dtatrc(self%m_size, 2, self%m_n_vars)) ! domrea.f90:216; TO DO: which shape?
        self%m_values_dtatrc(:, :, :) = huge(self%m_values_dtatrc(1, 1, 1)) ! domrea.f90:216
        allocate(self%m_values(self%m_size, self%m_n_vars)) ! domrea.f90:216
        self%m_values(:, :) = huge(self%m_values(1, 1)) ! domrea.f90:216

        self%m_alpha = alpha
        self%m_reduction_value_t = reduction_value_t
        self%m_length = length

    end subroutine init_members



    ! TO DO: check if it is true that the constructor has to be always overloaded
    ! TO DO: final version of the constructor should receive everything from a single namelist
    ! 'bc_name' is used just to avoid system used symbol 'name'

    !> Default constructor

    !> Calls bc default constructor and target constructor.
    type(sponge) function sponge_default(files_namelist, bc_name, n_vars, vars, var_names_idx, alpha, reduction_value_t, length)

        character(len=22), intent(in) :: files_namelist
        character(len=3) :: bc_name
        integer, intent(in) :: n_vars
        character(len=27), intent(in) :: vars
        integer(4), dimension(n_vars), intent(in) :: var_names_idx
        double precision, intent(in) :: alpha
        double precision, intent(in) :: reduction_value_t
        double precision, intent(in) :: length

        ! parent class constructor
        sponge_default%bc = bc(files_namelist)

        call sponge_default%init_members(bc_name, n_vars, vars, var_names_idx, alpha, reduction_value_t, length)

        write(*, *) 'INFO: successfully called sponge default constructor'

    end function sponge_default



    ! TO DO: check if it is true that the constructor has to be always overloaded
    ! TO DO: final version of the constructor should receive everything from a single namelist
    ! 'bc_name' is used just to avoid system used symbol 'name'

    !> Periodic constructor

    !> Calls bc periodic constructor and target constructor.
    type(sponge) function sponge_year(files_namelist, bc_name, n_vars, vars, var_names_idx, alpha, reduction_value_t, length, &
            start_time_string, end_time_string)

        character(len=27), intent(in) :: files_namelist
        character(len=3) :: bc_name
        integer, intent(in) :: n_vars
        character(len=27), intent(in) :: vars
        integer(4), dimension(n_vars), intent(in) :: var_names_idx
        double precision, intent(in) :: alpha
        double precision, intent(in) :: reduction_value_t
        double precision, intent(in) :: length
        character(len=17), intent(in) :: start_time_string
        character(len=17), intent(in) :: end_time_string

        ! parent class constructor
        sponge_year%bc = bc(files_namelist, start_time_string, end_time_string)

        call sponge_year%init_members(bc_name, n_vars, vars, var_names_idx, alpha, reduction_value_t, length)

        write(*, *) 'INFO: successfully called sponge year constructor'

    end function sponge_year



    !> Global size getter
    integer(4) function get_global_size(self)
        class(sponge), intent(in) :: self
        get_global_size = self%m_global_size
    end function get_global_size



    !> Overridden from bc
    subroutine load(self, idx)

        class(sponge), intent(inout) :: self
        integer, intent(in) :: idx
        integer :: i, j

        do i = 1, self%m_n_vars
            call readnc_double_1d(self%get_file_by_index(idx), self%m_var_names_data(i), self%m_global_size, self%m_aux)
            do j = 1, self%m_size
                self%m_values_dtatrc(j, 2, i) = self%m_aux(self%m_ridxt(1, j))
            enddo
        enddo

    end subroutine load



    !> Overridden from bc
    subroutine swap(self)

        class(sponge), intent(inout) :: self
        integer :: i, j

        do i = 1, self%m_n_vars
            do j = 1, self%m_size
                self%m_values_dtatrc(j, 1, i) = self%m_values_dtatrc(j, 2, i)
            enddo
        enddo

    end subroutine swap



    !> Overridden from bc
    subroutine actualize(self, weight)

        class(sponge), intent(inout) :: self
        double precision, intent(in) :: weight
        integer :: i, j

        do i = 1, self%m_n_vars
            do j = 1, self%m_size
                self%m_values(j, i) = (1.0 - weight) * self%m_values_dtatrc(j, 1, i) + weight * self%m_values_dtatrc(j, 2, i)
            enddo
        enddo

    end subroutine actualize



    !> Overridden from bc

    !> Actually it does not do anything, since the chosen policy is to use this class always in combination with a nudging.
    subroutine apply(self, e3t, n_tracers, trb, tra)

        use modul_param, only: jpk, jpj, jpi

        implicit none

        ! TO DO: to be removed. Find a way to enable both testing and production code.
        ! integer, parameter :: jpk = 70
        ! integer, parameter :: jpj = 65
        ! integer, parameter :: jpi = 182

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
        ! integer, parameter :: jpk = 70
        ! integer, parameter :: jpj = 65
        ! integer, parameter :: jpi = 182

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
                    idx_i = self%m_ridxt(4, j)
                    idx_j = self%m_ridxt(3, j)
                    idx_k = self%m_ridxt(2, j)
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
    subroutine apply_phys(self, lat, sponge_t, sponge_vel)

        use modul_param, only: jpk, jpj, jpi

        implicit none

        ! TO DO: to be removed. Find a way to enable both testing and production code.
        ! integer, parameter :: jpk = 70
        ! integer, parameter :: jpj = 65
        ! integer, parameter :: jpi = 182

        class(sponge), intent(inout) :: self
        double precision, dimension(jpj, jpi), intent(in) :: lat ! glamt
        double precision, dimension(jpj, jpi), intent(out) :: sponge_t ! spongeT
        double precision, dimension(jpk, jpj, jpi), intent(out) :: sponge_vel ! spongeVel

        integer :: i, j
        double precision :: reduction_value_vel ! reduction_value (only for sponge_vel)

        do i = 1, jpi
            do j = 1, jpj
                if (lat(j, i) < self%m_length) then
                    sponge_t(j, i) = self%m_reduction_value_t
                    reduction_value_vel = exp( -self%m_alpha * ((lat(j, i) - self%m_length)**2) )
                    sponge_vel(:, j, i) = reduction_value_vel
                endif
            enddo
        enddo

    end subroutine apply_phys



    !> Destructor
    subroutine sponge_destructor(self)

        class(sponge), intent(inout) :: self

        if (allocated(self%m_global_idxt)) then
            deallocate(self%m_global_idxt)
            write(*, *) 'INFO: m_global_idxt deallocated'
        endif

        if (allocated(self%m_var_names)) then
            deallocate(self%m_var_names)
            write(*, *) 'INFO: m_var_names deallocated'
        endif

        if (allocated(self%m_var_names_idxt)) then
            deallocate(self%m_var_names_idxt)
            write(*, *) 'INFO: m_var_names_idxt deallocated'
        endif

        if (allocated(self%m_var_names_data)) then
            deallocate(self%m_var_names_data)
            write(*, *) 'INFO: m_var_names_data deallocated'
        endif

        if (allocated(self%m_var_names_idx)) then
            deallocate(self%m_var_names_idx)
            write(*, *) 'INFO: m_var_names_idx deallocated'
        endif

        if (allocated(self%m_ridxt)) then
            deallocate(self%m_ridxt)
            write(*, *) 'INFO: m_ridxt deallocated'
        endif

        if (allocated(self%m_aux)) then
            deallocate(self%m_aux)
            write(*, *) 'INFO: m_aux deallocated'
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

    end subroutine sponge_destructor



end module sponge_mod
