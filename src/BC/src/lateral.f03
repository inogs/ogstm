module lateral_mod

    use boundary_mod

    implicit none

    private

    type, extends(boundary) :: lateral

        char(len=3) :: m_name ! ex: 'gib'; TO DO: should it belong to base class?

        char(len=15) :: m_bounmask ! 15 chars in order to handle names like 'bounmask_GIB.nc'

        integer :: m_n_vars ! BC_mem.f90:94
        char(len=3), allocatable, dimension(:) :: m_var_names ! domrea.F90:161-167
        char(len=5), allocatable, dimension(:) :: m_names_bounmask ! domrea.F90:172
        char(len=12), allocatable, dimension(:) :: m_names_idxt ! domrea.F90:204
        char(len=7), allocatable, dimension(:) :: m_names_data ! bc_gib.f90:113

        integer(4) :: m_global_size ! BC_mem.f90:20
        integer(4), allocatable, dimension(:) :: m_global_idxt ! BC_mem.f90:26
        integer(4) :: m_size ! BC_mem.f90:21

        ! TO DO other members missing

    contains

        private

        procedure :: set_global_size ! BC_mem.f90:70
        procedure :: set_global_idxt ! (call to readnc_int_1d(), domrea.f90:204-207,226)
        procedure :: count_in_subdomain ! (domrea.f90:209-211,228)
        procedure :: reindex ! (domrea.f90:218,233, which should be moved away with the 3D index)

        ! - put here the procedures to read file 'bounmask.nc' (domrea.F90:169-180)
        !   'bounmask.nc' anyway should not be needed, and all the parameters should be initialized inside the class

        ! delegated constructor
        procedure :: init_members ! (memory allocation also in domrea.f90:215-216,230:231)

        public

        final :: destructor

    end type lateral

    interface lateral
        module procedure lateral_default
        module procedure lateral_yearly
    end interface lateral

    public :: lateral

contains

    ! just a wrapper of 'getDimension' (BC_mem.f90:70)
    subroutine set_global_size(self)
        class(lateral), intent(in) :: self
        ! TO DO: 'get_file_by_index(integer index)' method to be implemented;
        !        remember to put a check: index must not excede time index
        !        (but it can excede file index...).
        call getDimension(self%m_data%get_file_by_index(1), self%m_names_idxt(1), self%m_global_size)
    end subroutine set_global_size

    ! just a wrapper of 'readnc_int_1d' (domrea.f90:204-207)
    subroutine set_global_idxt(self)
        class(lateral), intent(in) :: self
        call readnc_int_1d(self%m_data%get_file_by_index(1), self%m_names_idxt(1), self%m_global_size, self%m_global_idxt)
    end subroutine set_global_idxt

    subroutine init_members(self, name, bounmask, n_vars, vars)

        class(lateral), intent(in) :: self
        char(len=3) :: name
        char(len=15), intent(in) :: bounmask
        integer, intent(in) :: n_vars
        char(len=40), intent(in) :: vars ! max 10 vars, in the form: 'O2o N1p N3n N5s O3c O3h N6r ...'
        integer :: i

        self%m_name = name

        self%m_bounmask = bounmask

        self%m_n_vars = n_vars

        allocate(self%m_var_names(self%m_n_vars))
        allocate(self%m_names_bounmask(self%m_n_vars))
        allocate(self%m_names_idxt(self%m_n_vars))
        allocate(self%m_names_data(self%m_n_vars))

        do i = 1, self%m_n_vars
            integer :: start_idx, end_idx
            end_idx = 4*i - 1
            start_idx = end_idx - 2
            self%m_var_names(i) = vars(start_idx:end_idx)
            self%m_names_bounmask(i) = 're'//self%m_var_names(i)
            self%m_names_idxt(i) = self%m_name//'_idxt_'//self%m_var_names(i)
            slef%m_names_data(i) = self%m_name//'_'//self%m_var_names(i)
        enddo

        call self%set_global_size()

        allocate(self%m_global_idxt(self%m_global_size)) ! BC_mem.f90:111
        call self%set_global_idxt()

        ! TO DO: other members

    end subroutine init_members

    ! TO DO: check if it is true that the constructor has to be always overloaded
    ! TO DO: final version of the constructor should receive everything from a single namelist
    type(lateral) function lateral_default(files_namelist, name, bounmask, n_vars, vars)

        char(len=24), intent(in) :: files_namelist
        char(len=3) :: name
        char(len=15), intent(in) :: bounmask
        integer, intent(in) :: n_vars
        char(len=40), intent(in) :: vars

        lateral_default%boundary = boundary(files_namelist)

        call lateral_default%init_members(name, bounmask, n_vars, vars)

        ! ...

    end function lateral_default

    ! TO DO: check if it is true that the constructor has to be always overloaded
    ! TO DO: final version of the constructor should receive everything from a single namelist
    type(lateral) function lateral_yearly(files_namelist, name, bounmask, n_vars, vars, start_time_string, end_time_string)

        char(len=24), intent(in) :: files_namelist
        char(len=3) :: name
        char(len=15), intent(in) :: bounmask
        integer, intent(in) :: n_vars
        char(len=40), intent(in) :: vars
        char(len=17), intent(in) :: start_time_string
        char(len=17), intent(in) :: end_time_string

        lateral_yearly%boundary = boundary(files_namelist, start_time_string, end_time_string)

        call lateral_yearly%init_members(name, bounmask, n_vars, vars)

        ! ...

    end function lateral_yearly

    subroutine count_in_subdomain(self)
        class(lateral), intent(in) :: self
        ! ...
    end subroutine count_in_subdomain

    subroutine reindex(self)
        class(lateral), intent(in) :: self
        ! ...
    end subroutine reindex
    
    subroutine destructor(self)
        class(lateral), intent(in) :: self
    end subroutine

end module lateral_mod
