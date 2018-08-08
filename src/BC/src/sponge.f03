module sponge_mod

    use bc_mod
    use BC_mem

    implicit none

    private

    type, extends(bc) :: sponge

        char(len=3) :: m_name ! ex: 'gib'

        ! no more needed, bounmask should belong to nudging decorator
        ! char(len=15) :: m_bounmask ! 15 chars in order to handle names like 'bounmask_GIB.nc'

        integer(4) :: m_global_size ! BC_mem.f90:20
        integer(4), allocatable, dimension(:) :: m_global_idxt ! BC_mem.f90:26
        integer(4) :: m_size ! BC_mem.f90:21

        integer :: m_n_vars ! BC_mem.f90:94
        char(len=3), allocatable, dimension(:) :: m_var_names ! domrea.F90:161-167

        ! no more needed, bounmask should belong to nudging decorator
        ! char(len=5), allocatable, dimension(:) :: m_var_names_bounmask ! domrea.F90:172

        char(len=12), allocatable, dimension(:) :: m_var_names_idxt ! domrea.F90:204
        char(len=7), allocatable, dimension(:) :: m_var_names_data ! bc_gib.f90:113

        ! TO DO other members missing

    contains

        procedure :: set_global_size ! BC_mem.f90:70
        procedure :: set_global_idxt ! (call to readnc_int_1d(), domrea.f90:204-207,226)
        procedure :: count_in_subdomain ! (domrea.f90:209-211,228)
        procedure :: reindex ! (domrea.f90:218,233, which should be moved away with the 3D index)

        ! no more needed, bounmask should belong to nudging decorator
        ! put here the procedures to read file 'bounmask.nc' (domrea.F90:169-180)
        ! 'bounmask.nc' anyway should not be needed, and all the parameters should be initialized inside the class

        ! delegated constructor
        procedure :: init_members ! (memory allocation also in domrea.f90:215-216,230:231)

        final :: destructor

    end type sponge

    interface sponge
        module procedure sponge_default
        module procedure sponge_yearly
    end interface sponge

    public :: sponge

contains

    ! just a wrapper of 'getDimension' (BC_mem.f90:70)
    subroutine set_global_size(self)
        class(sponge), intent(in) :: self
        call getDimension(self%get_file_by_index(1), self%m_var_names_idxt(1), self%m_global_size)
    end subroutine set_global_size

    ! just a wrapper of 'readnc_int_1d' (domrea.f90:204-207)
    subroutine set_global_idxt(self)
        class(sponge), intent(in) :: self
        call readnc_int_1d(selfget_file_by_index(1), self%m_var_names_idxt(1), self%m_global_size, self%m_global_idxt)
    end subroutine set_global_idxt

    subroutine count_in_subdomain(self)
        class(sponge), intent(in) :: self
        self%m_size = COUNT_InSubDomain_GIB(self%m_global_size, self%m_global_idxt)
    end subroutine count_in_subdomain

    ! subroutine init_members(self, bc_name, bounmask, n_vars, vars)
    ! bc_name just to avoid system used symbol 'name'
    subroutine init_members(self, bc_name, n_vars, vars)

        class(sponge), intent(in) :: self
        char(len=3) :: bc_name
        ! char(len=15), intent(in) :: bounmask
        integer, intent(in) :: n_vars
        char(len=40), intent(in) :: vars ! max 10 vars, in the form: 'O2o N1p N3n N5s O3c O3h N6r ...'
        integer :: i

        self%m_name = bc_name

        ! self%m_bounmask = bounmask

        self%m_n_vars = n_vars

        allocate(self%m_var_names(self%m_n_vars))
        ! allocate(self%m_names_bounmask(self%m_n_vars))
        allocate(self%m_var_names_idxt(self%m_n_vars))
        allocate(self%m_var_names_data(self%m_n_vars))

        do i = 1, self%m_n_vars
            integer :: start_idx, end_idx
            end_idx = 4*i - 1
            start_idx = end_idx - 2
            self%m_var_names(i) = vars(start_idx:end_idx)
            ! self%m_names_bounmask(i) = 're'//self%m_var_names(i)
            self%m_var_names_idxt(i) = self%m_name//'_idxt_'//self%m_var_names(i)
            slef%m_var_names_data(i) = self%m_name//'_'//self%m_var_names(i)
        enddo

        call self%set_global_size()

        allocate(self%m_global_idxt(self%m_global_size)) ! BC_mem.f90:111

        call self%set_global_idxt()

        call self%count_in_subdomain()

        ! TO DO: other members

    end subroutine init_members

    ! TO DO: check if it is true that the constructor has to be always overloaded
    ! TO DO: final version of the constructor should receive everything from a single namelist
    ! type(sponge) function sponge_default(files_namelist, name, bounmask, n_vars, vars)
    ! bc_name just to avoid system used symbol 'name'
    type(sponge) function sponge_default(files_namelist, bc_name, n_vars, vars)

        char(len=24), intent(in) :: files_namelist
        char(len=3) :: bc_name
        ! char(len=15), intent(in) :: bounmask
        integer, intent(in) :: n_vars
        char(len=40), intent(in) :: vars

        sponge_default%bc = bc(files_namelist)

        ! call sponge_default%init_members(name, bounmask, n_vars, vars)
        call sponge_default%init_members(bc_name, n_vars, vars)

        ! ...

    end function sponge_default

    ! TO DO: check if it is true that the constructor has to be always overloaded
    ! TO DO: final version of the constructor should receive everything from a single namelist
    ! bc_name just to avoid system used symbol 'name'
    type(sponge) function sponge_yearly(files_namelist, bc_name, bounmask, n_vars, vars, start_time_string, end_time_string)

        char(len=24), intent(in) :: files_namelist
        char(len=3) :: bc_name
        char(len=15), intent(in) :: bounmask
        integer, intent(in) :: n_vars
        char(len=40), intent(in) :: vars
        char(len=17), intent(in) :: start_time_string
        char(len=17), intent(in) :: end_time_string

        sponge_yearly%bc = bc(files_namelist, start_time_string, end_time_string)

        call sponge_yearly%init_members(bc_name, bounmask, n_vars, vars)

        ! ...

    end function sponge_yearly

    subroutine reindex(self)
        class(sponge), intent(in) :: self
        ! ...
    end subroutine reindex
    
    subroutine destructor(self)
        class(sponge), intent(in) :: self
    end subroutine

end module sponge_mod
