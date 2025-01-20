!> Base class from which every boundary types must inherit.
module bc_mod

    use bc_data_mod

    implicit none

    private

    type bc
        ! Eventually this should became an array of m_bc_data, one for each variable
        type(bc_data), pointer :: m_bc_data => null()
    contains
        procedure :: const_data
        procedure :: get_file_by_index
        procedure :: get_prev_idx
        procedure :: get_next_idx
        procedure :: set_current_interval
        procedure :: get_interpolation_factor
        procedure :: load
        procedure :: swap
        procedure :: actualize
        procedure :: update_device
        procedure :: apply
        procedure :: apply_dirichlet
        procedure :: apply_nudging
        procedure :: apply_phys
        procedure :: fix_diagnostic_vars ! overridden only by hard open class
        procedure :: bc_destructor
    end type bc

    interface bc
        module procedure bc_empty
        module procedure bc_default
        module procedure bc_year
    end interface bc

    public :: bc

contains



    !> bc empty constructor
    type(bc) function bc_empty()

        allocate(bc_empty%m_bc_data)
        bc_empty%m_bc_data = bc_data()

        ! write(*, *) 'INFO: successfully called bc empty constructor'

    end function bc_empty



    !> bc default constructor.

    !> Calls bc_data default contructor.
    type(bc) function bc_default(filenames_list)

        character(len=22), intent(in) :: filenames_list

        allocate(bc_default%m_bc_data)
        bc_default%m_bc_data = bc_data(filenames_list)

        ! write(*, *) 'INFO: successfully called bc default constructor'

    end function bc_default



    !> bc periodic constructor

    !> Calls bc_data periodic constructor
    type(bc) function bc_year(filenames_list, start_time_string, end_time_string)

        character(len=22), intent(in) :: filenames_list
        character(len=17), intent(in) :: start_time_string
        character(len=17), intent(in) :: end_time_string

        allocate(bc_year%m_bc_data)
        bc_year%m_bc_data = bc_data(filenames_list, start_time_string, end_time_string)

        ! write(*, *) 'INFO: successfully called bc empty constructor'

    end function bc_year



    !> Constant data getter

    !> Boolean getter to check whether data are constant or not
    logical function const_data(self)
        class(bc), intent(in) :: self
        const_data = self%m_bc_data%const()
    end function const_data



    !> Getter for the data file given the time index
    character(len=100) function get_file_by_index(self, idx)

        class(bc), intent(in) :: self
        integer, intent(in) :: idx

        get_file_by_index = self%m_bc_data%get_file_by_index(idx)

    end function get_file_by_index



    !> Getter for the left extreme of the current time interval
    integer function get_prev_idx(self)
        class(bc), intent(in) :: self
        get_prev_idx = self%m_bc_data%get_prev_idx()
    end function get_prev_idx



    !> Getter for the right extreme of the current time interval
    integer function get_next_idx(self)
        class(bc), intent(in) :: self
        get_next_idx = self%m_bc_data%get_next_idx()
    end function get_next_idx



    !> Setter of the current time interval

    !> Just sets the current time interval (prev and next indexes) according to current_time_string,
    !! without computing the interpolation factor.
    subroutine set_current_interval(self, current_time_string, new_data)

        class(bc), intent(inout) :: self
        character(len=17), intent(in) :: current_time_string
        logical, optional, intent(out) :: new_data

        call self%m_bc_data%set_current_interval(current_time_string)
        if (present(new_data)) then
            new_data = self%m_bc_data%new_interval()
        endif

    end subroutine set_current_interval



    !> Linear interpolation factor computation

    !> Computes and returns the linear interpolation factor,
    !! keeping track of the current time interval (prev and next indexes).
    double precision function get_interpolation_factor(self, current_time_string, new_data)

        class(bc), intent(inout) :: self
        character(len=17), intent(in) :: current_time_string
        logical, optional, intent(out) :: new_data

        get_interpolation_factor = self%m_bc_data%get_interpolation_factor(current_time_string)
        if (present(new_data)) then
            new_data = self%m_bc_data%new_interval()
        endif

    end function get_interpolation_factor



    !> Loads in memory the data of the extreme of a time interval.
    subroutine load(self, idx)
        class(bc), intent(inout) :: self
        integer, intent(in) :: idx
        write(*, *) 'WARN: base class does not implement method load'
    end subroutine load



    !> Swaps in memory the data of two extremes.
    subroutine swap(self)
        class(bc), intent(inout) :: self
        write(*, *) 'WARN: base class does not implement method swap'
    end subroutine swap



    !> Sets the right values according to the interpolation weight

    !> Requires a weight to correct the loaded data,
    !! even if the weight of the linear interpolation can be computed by the class itself
    !! through its bc_data object.
    !! In this way the method is more flexible
    !! when the linear interpolation option is not enabled in the namelist file,
    !! and the weight needs to be set always to 1.
    subroutine actualize(self, weight)
        class(bc), intent(inout) :: self
        double precision, intent(in) :: weight
        write(*, *) 'WARN: base class does not implement method actualize'
    end subroutine actualize

    subroutine update_device(self)
        class(bc), intent(inout) :: self
        write(*, *) 'WARN: base class does not implement method update_device'
    end subroutine update_device

    !> Applies the tracer values at the boundaries to the final tracers matrix.

    !> Called to set the final values of the tracer fields near the boundaries,
    !! according to the type of boundary.
    !! Model output fields are adjusted according to the boundary scheme.
    subroutine apply(self, e3t, n_tracers, trb, tra)

        use modul_param, only: jpk, jpj, jpi

        implicit none

        ! TO DO: to be removed. Find a way to enable both testing and production code.
        ! integer, parameter :: jpk = 125
        ! integer, parameter :: jpj = 380
        ! integer, parameter :: jpi = 1085

        class(bc), intent(inout) :: self
        double precision, dimension(jpk, jpj, jpi), intent(in) :: e3t
        integer, intent(in) :: n_tracers
        double precision, dimension(jpk, jpj, jpi, n_tracers), intent(in) :: trb
        double precision, dimension(jpk, jpj, jpi, n_tracers), intent(inout) :: tra

        write(*, *) 'WARN: base class does not implement method apply'

    end subroutine apply

    subroutine apply_dirichlet(self)

        use modul_param, only: jpk, jpj, jpi
        use myalloc

        implicit none

        ! TO DO: to be removed. Find a way to enable both testing and production code.

        class(bc), intent(inout) :: self

        write(*, *) 'WARN bc.f90: base class does not implement method apply_dirichlet'

    end subroutine apply_dirichlet

    !> Used to simplify the apply method if a nudging scheme has to be included in the fields correction.
    subroutine apply_nudging(self, e3t, n_tracers, rst_tracers, trb, tra)

        use modul_param, only: jpk, jpj, jpi

        implicit none

        ! TO DO: to be removed. Find a way to enable both testing and production code.
        ! integer, parameter :: jpk = 125
        ! integer, parameter :: jpj = 380
        ! integer, parameter :: jpi = 1085

        class(bc), intent(inout) :: self
        double precision, dimension(jpk, jpj, jpi), intent(in) :: e3t
        integer, intent(in) :: n_tracers
        double precision, dimension(jpk, jpj, jpi, n_tracers), intent(in) :: rst_tracers
        double precision, dimension(jpk, jpj, jpi, n_tracers), intent(in) :: trb
        double precision, dimension(jpk, jpj, jpi, n_tracers), intent(inout) :: tra

        write(*, *) 'WARN: base class does not implement method apply_nudging'

    end subroutine apply_nudging



    !> Updates the values of the OGCM at the boundaries.

    !> Used to adjust the values of the velocity fields at the boundaries.
    !! Velocity fields are imposed to the model,
    !! usually according to the output of the NEMO Ocean Model or of an equivalent OGCM.
    !! Anyway they are not computed directly by OGSTM.
    !! Nevertheless with some specific boundary conditions schemes these outputs still need to be modified at the boundaries.
    !! This is the case for example of the sponge boundary,
    !! in which the velocity fields are forcefully set to zero at the boundary of OGSTM,
    !! in order for it to be closed.
    !! In such cases, this method is necessary both
    !! to nullify the velocity field components and
    !! to smooth the actual values of the OGCM forcing fields.
    subroutine apply_phys(self, lon, sponge_t, sponge_vel)

        use modul_param, only: jpk, jpj, jpi

        implicit none

        ! TO DO: to be removed. Find a way to enable both testing and production code.
        ! integer, parameter :: jpk = 125
        ! integer, parameter :: jpj = 380
        ! integer, parameter :: jpi = 1085

        class(bc), intent(inout) :: self
        double precision, dimension(jpj, jpi), intent(in) :: lon
        double precision, dimension(jpj, jpi), intent(out) :: sponge_t
        double precision, dimension(jpk, jpj, jpi), intent(out) :: sponge_vel

        write(*, *) 'WARN: base class does not implement method apply_phys'

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

        class(bc), intent(inout) :: self

        ! write(*, *) 'WARN: base class does not implement method fix_diagnostic_vars. Only hard open overrides it'

    end subroutine fix_diagnostic_vars



    !> Destructor
    subroutine bc_destructor(self)

        class(bc), intent(inout) :: self

        ! First call bc_data destructor
        call self%m_bc_data%bc_data_destructor()
        
        deallocate(self%m_bc_data)
        ! write(*, *) 'INFO: m_bc_data deallocated'
        nullify(self%m_bc_data)
        ! write(*, *) 'INFO: m_bc_data deassociated'

    end subroutine bc_destructor



end module bc_mod
