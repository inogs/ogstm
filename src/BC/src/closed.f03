module closed_mod

    use bc_mod
    use bc_aux_mod

    implicit none

    private

    type, extends(bc) :: closed

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
        procedure :: closed_destructor

    end type closed

    interface closed
        module procedure closed_default
        module procedure closed_year
    end interface closed

    public :: closed

contains



    ! just a wrapper of 'getDimension' (BC_mem.f90:70)
    subroutine set_global_size(self)
        class(closed), intent(inout) :: self
        call getDimension(self%get_file_by_index(1), self%m_var_names_idxt(1), self%m_global_size)
    end subroutine set_global_size



    ! just a wrapper of 'readnc_int_1d' (domrea.f90:204-207)
    subroutine set_global_idxt(self)
        class(closed), intent(inout) :: self
        allocate(self%m_global_idxt(self%m_global_size)) ! BC_mem.f90:111
        self%m_global_idxt(:) = huge(self%m_global_idxt(1))
        call readnc_int_1d(self%get_file_by_index(1), self%m_var_names_idxt(1), self%m_global_size, self%m_global_idxt)
    end subroutine set_global_idxt



    ! just a wrapper of 'COUNT_InSubDomain_GIB' (domrea.f90:209)
    subroutine set_size(self)
        class(closed), intent(inout) :: self
        self%m_size = COUNT_InSubDomain_3d(self%m_global_size, self%m_global_idxt)
    end subroutine set_size



    ! just a wrapper of 'GIBRE_Indexing' (domrea.f90:218)
    subroutine reindex(self)
        class(closed), intent(inout) :: self
        call RE_Indexing_3d(self%m_global_size, self%m_global_idxt, self%m_size, self%m_ridxt)
    end subroutine reindex



    ! subroutine init_members(self, bc_name, bounmask, n_vars, vars)
    ! 'bc_name' is used just to avoid system used symbol 'name'
    subroutine init_members(self, bc_name, n_vars, vars, var_names_idx)

        class(closed), intent(inout) :: self
        character(len=3) :: bc_name
        integer, intent(in) :: n_vars
        character(len=27), intent(in) :: vars ! 'O2o N1p N3n N5s O3c O3h N6r'; TO DO: more flexible
        integer(4), dimension(n_vars), intent(in) :: var_names_idx
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

    end subroutine init_members



    ! TO DO: check if it is true that the constructor has to be always overloaded
    ! TO DO: final version of the constructor should receive everything from a single namelist
    ! 'bc_name' is used just to avoid system used symbol 'name'
    type(closed) function closed_default(files_namelist, bc_name, n_vars, vars, var_names_idx)

        character(len=22), intent(in) :: files_namelist
        character(len=3) :: bc_name
        integer, intent(in) :: n_vars
        character(len=27), intent(in) :: vars
        integer(4), dimension(n_vars), intent(in) :: var_names_idx

        ! parent class constructor
        closed_default%bc = bc(files_namelist)

        call closed_default%init_members(bc_name, n_vars, vars, var_names_idx)

        write(*, *) 'INFO: successfully called closed default constructor'

    end function closed_default



    ! TO DO: check if it is true that the constructor has to be always overloaded
    ! TO DO: final version of the constructor should receive everything from a single namelist
    ! 'bc_name' is used just to avoid system used symbol 'name'
    type(closed) function closed_year(files_namelist, bc_name, n_vars, vars, var_names_idx, start_time_string, end_time_string)

        character(len=27), intent(in) :: files_namelist
        character(len=3) :: bc_name
        integer, intent(in) :: n_vars
        character(len=27), intent(in) :: vars
        integer(4), dimension(n_vars), intent(in) :: var_names_idx
        character(len=17), intent(in) :: start_time_string
        character(len=17), intent(in) :: end_time_string

        ! parent class constructor
        closed_year%bc = bc(files_namelist, start_time_string, end_time_string)

        call closed_year%init_members(bc_name, n_vars, vars, var_names_idx)

        write(*, *) 'INFO: successfully called closed year constructor'

    end function closed_year



    integer(4) function get_global_size(self)
        class(closed), intent(in) :: self
        get_global_size = self%m_global_size
    end function get_global_size



    subroutine load(self, idx)

        class(closed), intent(inout) :: self
        integer, intent(in) :: idx
        integer :: i, j

        do i = 1, self%m_n_vars
            call readnc_double_1d(self%get_file_by_index(idx), self%m_var_names_data(i), self%m_global_size, self%m_aux)
            do j = 1, self%m_size
                self%m_values_dtatrc(j, 2, i) = self%m_aux(self%m_ridxt(1, j))
            enddo
        enddo

    end subroutine load



    subroutine swap(self)

        class(closed), intent(inout) :: self
        integer :: i, j

        do i = 1, self%m_n_vars
            do j = 1, self%m_size
                self%m_values_dtatrc(j, 1, i) = self%m_values_dtatrc(j, 2, i)
            enddo
        enddo

    end subroutine swap



    subroutine actualize(self, weight)

        class(closed), intent(inout) :: self
        double precision, intent(in) :: weight
        integer :: i, j

        do i = 1, self%m_n_vars
            do j = 1, self%m_size
                self%m_values(j, i) = (1.0 - weight) * self%m_values_dtatrc(j, 1, i) + weight * self%m_values_dtatrc(j, 2, i)
            enddo
        enddo

    end subroutine actualize



    subroutine apply(self, e3t, n_tracers, trb, tra)

        use modul_param, only: jpk, jpj, jpi

        implicit none

        ! TO DO: to be removed. Find a way to enable both testing and production code.
        ! integer, parameter :: jpk = 70
        ! integer, parameter :: jpj = 65
        ! integer, parameter :: jpi = 182

        class(closed), intent(inout) :: self
        double precision, dimension(jpk, jpj, jpi), intent(in) :: e3t
        integer, intent(in) :: n_tracers
        double precision, dimension(jpk, jpj, jpi, n_tracers), intent(in) :: trb
        double precision, dimension(jpk, jpj, jpi, n_tracers), intent(inout) :: tra

        write(*, *) 'WARN: closed class does not implement this method'
        write(*, *) 'WARN: attempt to apply boundary conditions on closed boundary without nudging'

    end subroutine apply



    subroutine apply_nudging(self, e3t, n_tracers, rst_tracers, trb, tra)

        use modul_param, only: jpk, jpj, jpi

        implicit none

        ! TO DO: to be removed. Find a way to enable both testing and production code.
        ! integer, parameter :: jpk = 70
        ! integer, parameter :: jpj = 65
        ! integer, parameter :: jpi = 182

        class(closed), intent(inout) :: self
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



    ! TO DO: provide more flexibility both
    ! - in the choice of the 'closed' function
    ! - in the choice of the interested area
    ! So far, only a rectangular area in the west boundary has been implemented
    ! WARNING: closed_t and closed_vel will be overwritten every time the method is called!
    subroutine apply_phys(self, lat, sponge_t, sponge_vel)

        use modul_param, only: jpk, jpj, jpi

        implicit none

        ! TO DO: to be removed. Find a way to enable both testing and production code.
        ! integer, parameter :: jpk = 70
        ! integer, parameter :: jpj = 65
        ! integer, parameter :: jpi = 182

        class(closed), intent(inout) :: self
        double precision, dimension(jpj, jpi), intent(in) :: lat ! glamt
        double precision, dimension(jpj, jpi), intent(out) :: sponge_t ! spongeT
        double precision, dimension(jpk, jpj, jpi), intent(out) :: sponge_vel ! spongeVel

        write(*, *) 'WARN: sponge_t and sponge_val are left untouched by closed class'

    end subroutine apply_phys



    subroutine closed_destructor(self)

        class(closed), intent(inout) :: self

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

    end subroutine closed_destructor



end module closed_mod
