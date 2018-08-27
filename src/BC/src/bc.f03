module bc_mod

    use bc_data_mod

    implicit none

    private

    type bc
        ! Eventually this should became an array of m_bc_data, one for each variable
        type(bc_data), pointer :: m_bc_data => null()
    contains
        procedure :: get_file_by_index
        procedure :: load
        procedure :: swap
        procedure :: actualize
        procedure :: apply
        procedure :: apply_phys
        procedure :: bc_destructor
    end type bc

    interface bc
        module procedure bc_empty
        module procedure bc_default
        module procedure bc_year
    end interface bc

    public :: bc

contains



    type(bc) function bc_empty()

        allocate(bc_empty%m_bc_data)
        bc_empty%m_bc_data = bc_data()

    end function bc_empty



    type(bc) function bc_default(files_namelist)

        character(len=22), intent(in) :: files_namelist

        allocate(bc_default%m_bc_data)
        bc_default%m_bc_data = bc_data(files_namelist)

    end function bc_default



    type(bc) function bc_year(files_namelist, start_time_string, end_time_string)

        character(len=27), intent(in) :: files_namelist
        character(len=17), intent(in) :: start_time_string
        character(len=17), intent(in) :: end_time_string

        allocate(bc_year%m_bc_data)
        bc_year%m_bc_data = bc_data(files_namelist, start_time_string, end_time_string)

    end function bc_year



    character(len=24) function get_file_by_index(self, idx)

        class(bc), intent(in) :: self
        integer, intent(in) :: idx

        get_file_by_index = self%m_bc_data%get_file_by_index(idx)

    end function get_file_by_index



    subroutine load(self, idx)
        class(bc), intent(inout) :: self
        integer, intent(in) :: idx
        write(*, *) 'WARN: base class does not implement this method'
    end subroutine load



    subroutine swap(self)
        class(bc), intent(inout) :: self
        write(*, *) 'WARN: base class does not implement this method'
    end subroutine swap



    subroutine actualize(self, weight)
        class(bc), intent(inout) :: self
        double precision, intent(in) :: weight
        write(*, *) 'WARN: base class does not implement this method'
    end subroutine actualize



    subroutine apply(self, e3t, n_tracers, rst_tracers, trb, tra)

        use modul_param, only: jpk, jpj, jpi

        implicit none

        ! TO DO: to be removed. Find a way to enable both testing and production code.
        ! integer, parameter :: jpk = 70
        ! integer, parameter :: jpj = 65
        ! integer, parameter :: jpi = 182

        class(bc), intent(inout) :: self
        double precision, dimension(jpk, jpj, jpi), intent(in) :: e3t
        integer, intent(in) :: n_tracers
        double precision, dimension(jpk, jpj, jpi, n_tracers), intent(in) :: rst_tracers
        double precision, dimension(jpk, jpj, jpi, n_tracers), intent(in) :: trb
        double precision, dimension(jpk, jpj, jpi, n_tracers), intent(inout) :: tra

        write(*, *) 'WARN: base class does not implement this method'

    end subroutine apply



    subroutine apply_phys(self, lat, sponge_t, sponge_vel)

        use modul_param, only: jpk, jpj, jpi

        implicit none

        ! TO DO: to be removed. Find a way to enable both testing and production code.
        ! integer, parameter :: jpk = 70
        ! integer, parameter :: jpj = 65
        ! integer, parameter :: jpi = 182

        class(bc), intent(inout) :: self
        double precision, dimension(jpj, jpi), intent(in) :: lat
        double precision, dimension(jpj, jpi), intent(out) :: sponge_t
        double precision, dimension(jpk, jpj, jpi), intent(out) :: sponge_vel

        write(*, *) 'WARN: base class does not implement this method'

    end subroutine apply_phys



    subroutine bc_destructor(self)

        class(bc), intent(inout) :: self

        ! First call bc_data destructor
        call self%m_bc_data%bc_data_destructor()
        
        deallocate(self%m_bc_data)
        write(*, *) 'INFO: m_bc_data deallocated'
        nullify(self%m_bc_data)
        write(*, *) 'INFO: m_bc_data deassociated'

    end subroutine bc_destructor



end module bc_mod
