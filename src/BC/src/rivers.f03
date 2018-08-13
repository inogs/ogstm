module rivers_mod

    use bc_mod
    use bc_aux_mod

    implicit none

    private

    type, extends(bc) :: rivers

        character(len=3) :: m_name ! ex: 'riv'
        integer(4) :: m_global_size ! BC_mem.f90:20
        integer(4), allocatable, dimension(:) :: m_global_idxt ! BC_mem.f90:26; TO DO: find better name
        integer(4) :: m_size ! BC_mem.f90:21
        integer :: m_n_vars ! BC_mem.f90:95
        character(len=3), allocatable, dimension(:) :: m_var_names
        character(len=12) :: m_name_idxt ! domrea.f90:226; TO DO: find better name
        character(len=7), allocatable, dimension(:) :: m_var_names_data ! bc_tin.f90:116
        integer(4), allocatable, dimension(:, :) :: m_ridxt ! TO DO: find better name
        double precision, allocatable, dimension(:) :: m_aux
        double precision, allocatable, dimension(:, :, :) :: m_values_dtatrc ! TO DO: find better name
        double precision, allocatable, dimension(:, :) :: m_values

    contains

        ! delegated constructor - related procedures
        procedure :: set_global_size ! BC_mem.f90:71
        procedure :: set_global_idxt ! (call to readnc_int_1d(), domrea.f90:226)
        procedure :: set_size ! (domrea.f90:228)
        procedure :: reindex ! (domrea.f90:233, which should be moved away with the 3D index)
        ! delegated constructor
        procedure :: init_members ! (memory allocation also in domrea.f90:230:231)
        ! getters
        procedure :: get_global_size
        ! base class methods
        procedure :: load
        procedure :: swap
        procedure :: actualize
        ! destructor
        procedure :: rivers_destructor

    end type rivers

    interface rivers
        module procedure rivers_default
        module procedure rivers_year
    end interface rivers

    public :: rivers

contains



    ! just a wrapper of 'getDimension' (BC_mem.f90:71)
    subroutine set_global_size(self)
        class(rivers), intent(inout) :: self
        call getDimension(self%get_file_by_index(1), self%m_var_names_idxt(1), self%m_global_size)
    end subroutine set_global_size



    ! just a wrapper of 'readnc_int_1d' (domrea.f90:226)
    subroutine set_global_idxt(self)
        class(rivers), intent(inout) :: self
        allocate(self%m_global_idxt(self%m_global_size)) ! BC_mem.f90:135
        self%m_global_idxt(:) = huge(self%m_global_idxt(1))
        call readnc_int_1d(self%get_file_by_index(1), self%m_name_idxt, self%m_global_size, self%m_global_idxt)
    end subroutine set_global_idxt



    ! just a wrapper of 'COUNT_InSubDomain' (domrea.f90:228)
    subroutine set_size(self)
        class(rivers), intent(inout) :: self
        self%m_size = COUNT_InSubDomain(self%m_global_size, self%m_global_idxt)
    end subroutine set_size



    ! just a wrapper of 'RIVRE_Indexing' (domrea.f90:233)
    subroutine reindex(self)
        class(rivers), intent(inout) :: self
        call RE_Indexing(self%m_global_size, self%m_global_idxt, self%m_size, self%m_ridxt)
    end subroutine reindex



    ! 'bc_name' is used just to avoid system used symbol 'name'
    subroutine init_members(self, bc_name, n_vars, vars)

        class(rivers), intent(inout) :: self
        character(len=3) :: bc_name
        integer, intent(in) :: n_vars
        character(len=23), intent(in) :: vars ! 'N3n N1p N5s O3c O3h O2o'; TO DO: more flexible
        integer :: i, start_idx, end_idx

        self%m_name = bc_name
        self%m_n_vars = n_vars
        self%m_name_idxt = self%m_name//'_idxt'

        allocate(self%m_var_names(self%m_n_vars))
        allocate(self%m_var_names_data(self%m_n_vars))

        do i = 1, self%m_n_vars
            end_idx = 4*i - 1
            start_idx = end_idx - 2
            self%m_var_names(i) = vars(start_idx:end_idx)
            self%m_var_names_data(i) = self%m_name//'_'//self%m_var_names(i)
        enddo

        ! call delegated constructor - related procedures
        call self%set_global_size()
        call self%set_global_idxt()
        call self%set_size()

        allocate(self%m_ridxt(4, self%m_size)) ! domrea.f90:231
        self%m_ridxt(:, :) = huge(self%m_ridxt(1, 1)) ! domrea.f90:231
        call self%reindex() ! domrea.f90:233

        allocate(self%m_aux(self%m_global_size))
        self%m_aux(:) = huge(self%m_aux(1))
        allocate(self%m_values_dtatrc(2, self%m_size, self%m_n_vars)) ! domrea.f90:216; TO DO: which shape?
        self%m_values_dtatrc(:, :, :) = huge(self%m_values_dtatrc(1, 1, 1)) ! domrea.f90:231
        allocate(self%m_values(self%m_size, self%m_n_vars)) ! domrea.f90:231
        self%m_values(:, :) = huge(self%m_values(1, 1)) ! domrea.f90:231

    end subroutine init_members



    ! TO DO: check if it is true that the constructor has to be always overloaded
    ! TO DO: final version of the constructor should receive everything from a single namelist
    ! 'bc_name' is used just to avoid system used symbol 'name'
    type(rivers) function rivers_default(files_namelist, bc_name, n_vars, vars)

        character(len=22), intent(in) :: files_namelist
        character(len=3) :: bc_name
        integer, intent(in) :: n_vars
        character(len=23), intent(in) :: vars

        ! parent class constructor
        rivers_default%bc = bc(files_namelist)

        call rivers_default%init_members(bc_name, n_vars, vars)

    end function rivers_default



    ! TO DO: check if it is true that the constructor has to be always overloaded
    ! TO DO: final version of the constructor should receive everything from a single namelist
    ! 'bc_name' is used just to avoid system used symbol 'name'
    type(rivers) function rivers_year(files_namelist, bc_name, n_vars, vars, start_time_string, end_time_string)

        character(len=27), intent(in) :: files_namelist
        character(len=3) :: bc_name
        integer, intent(in) :: n_vars
        character(len=23), intent(in) :: vars
        character(len=17), intent(in) :: start_time_string
        character(len=17), intent(in) :: end_time_string

        ! parent class constructor
        rivers_year%bc = bc(files_namelist, start_time_string, end_time_string)

        call rivers_year%init_members(bc_name, n_vars, vars)

    end function rivers_year



    integer(4) function get_global_size(self)
        class(rivers), intent(in) :: self
        get_global_size = self%m_global_size
    end function get_global_size



    subroutine load(self, idx)

        class(rivers), intent(inout) :: self
        integer, intent(in) :: idx
        integer :: i, j

        do i = 1, self%m_n_vars
            call readnc_double_1d(self%get_file_by_index(idx), self%m_var_names_data(i), self%m_global_size, self%m_aux)
            do j = 1, self%m_size
                self%m_values_dtatrc(2, j, i) = self%m_aux(self%m_ridxt(1, j))
            enddo
        enddo

    end subroutine load



    subroutine swap(self)

        class(rivers), intent(inout) :: self
        integer :: i, j

        do i = 1, self%m_n_vars
            do j = 1, self%m_size
                self%m_values_dtatrc(1, j, i) = self%m_values_dtatrc(2, j, i)
            enddo
        enddo

    end subroutine swap



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



    subroutine rivers_destructor(self)

        class(rivers), intent(inout) :: self

        if (allocated(self%m_global_idxt)) then
            deallocate(self%m_global_idxt)
            write(*, *) 'INFO: m_global_idxt deallocated'
        endif

        if (allocated(self%m_var_names)) then
            deallocate(self%m_var_names)
            write(*, *) 'INFO: m_var_names deallocated'
        endif

        if (allocated(self%m_var_names_data)) then
            deallocate(self%m_var_names_data)
            write(*, *) 'INFO: m_var_names_data deallocated'
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

    end subroutine rivers_destructor



end module rivers_mod
