!> Class bc_set represents a set of boundary conditions

!> It is supposed to be used to store and operate on all the boundary conditions imposed to the model.
!! It features an array of pointers to bc class objects, containing all the boundary conditions.
module bc_set_mod

    use myalloc, only: lwp
    use bc_mod
    use rivers_mod
    use sponge_mod
    use hard_open_mod
    use nudging_mod
    use bc_handle_mod

    implicit none

    private

    type bc_container
        class(bc), pointer :: content => null()
    end type bc_container

    type bc_set
        integer :: m_n_bcs
        type(bc_container), allocatable, dimension(:) :: m_bcs
    contains
        procedure :: update
        procedure :: apply
        procedure :: apply_dirichlet
        procedure :: apply_phys
        procedure :: fix_diagnostic_vars
        procedure :: bc_set_destructor
    end type bc_set

    interface bc_set
        module procedure bc_set_default
    end interface bc_set

    public :: bc_set

contains



    subroutine bc_set_default(bcs,bcs_namelist)

        type(bc_set), intent(inout) :: bcs
        character(len=14), intent(in) :: bcs_namelist
        integer, parameter :: file_unit = 102 ! 100 for data files, 101 for boundary namelist files, 102 for global namelist
        character(len=54) :: bc_string
        integer :: i

        ! open file
        open(unit=file_unit, file=bcs_namelist)
        if (lwp) write(*, *) 'INFO: reading from file ', bcs_namelist

        ! get number of boundaries and allocate memory accordingly
        read(file_unit, *) bcs%m_n_bcs
        allocate(bcs%m_bcs(bcs%m_n_bcs))

        ! get info for each boundary and initialize boundaries accordingly
        do i = 1, bcs%m_n_bcs
            read(file_unit, *) bc_string
            if (lwp) write(*, *) 'INFO: initializing boundary with bc_string = ', bc_string
            bcs%m_bcs(i)%content => bc_init(bc_string)
        enddo

        ! close file
        close(unit=file_unit)

    end subroutine bc_set_default



    subroutine update(self, current_time_string)

        class(bc_set), intent(inout) :: self
        character(len=17), intent(in) :: current_time_string
        integer :: i

        do i = 1, self%m_n_bcs
            call bc_update(self%m_bcs(i)%content, current_time_string)
        enddo

    end subroutine update



    subroutine apply(self, e3t, trb, tra)

        use modul_param, only: jpk, jpj, jpi

        implicit none

        class(bc_set), intent(inout) :: self
        double precision, dimension(jpk, jpj, jpi), intent(in) :: e3t
        double precision, dimension(jpk, jpj, jpi, jptra), intent(in) :: trb
        double precision, dimension(jpk, jpj, jpi, jptra), intent(inout) :: tra
        integer :: i

        do i = 1, self%m_n_bcs
            call self%m_bcs(i)%content%apply(e3t, jptra, trb, tra)
        enddo

    end subroutine apply

    subroutine apply_dirichlet(self)

        use modul_param, only: jpk, jpj, jpi
        use myalloc
        implicit none
        integer i

        class(bc_set), intent(inout) :: self
        class(bc), pointer :: bc_ptr ! auxiliary pointer to guarantee polymorphism

        do i = 1,self%m_n_bcs
           bc_ptr => self%m_bcs(i)%content
           SELECT TYPE(bc_ptr)
               CLASS IS (hard_open)
                  call self%m_bcs(i)%content%apply_dirichlet()
           END SELECT
        enddo



    end subroutine apply_dirichlet

    subroutine apply_phys(self, lon, sponge_t, sponge_vel, internal_sponging)

        use modul_param, only: jpk, jpj, jpi

        implicit none

        class(bc_set), intent(inout) :: self
        double precision, dimension(jpj, jpi), intent(in) :: lon
        double precision, dimension(jpj, jpi), intent(out) :: sponge_t
        double precision, dimension(jpk, jpj, jpi), intent(out) :: sponge_vel
        logical, intent(in) :: internal_sponging ! internal_sponging is set from global namelist file
        integer :: i

        if (internal_sponging) then
            do i = 1, self%m_n_bcs
                call self%m_bcs(i)%content%apply_phys(lon, sponge_t, sponge_vel)
            enddo
        endif

    end subroutine apply_phys



    subroutine fix_diagnostic_vars(self, tra_dia, tra_dia_2d)

        use modul_param, only: jpk, jpj, jpi

        implicit none

        class(bc_set), intent(inout) :: self
        double precision, dimension(jptra_dia, jpk, jpj, jpi), intent(inout) :: tra_dia
        double precision, dimension(jptra_dia_2d, jpj, jpi), intent(inout) :: tra_dia_2d
        integer :: i

        do i = 1, self%m_n_bcs
            call self%m_bcs(i)%content%fix_diagnostic_vars(jptra_dia, tra_dia, jptra_dia_2d, tra_dia_2d)
        enddo

    end subroutine fix_diagnostic_vars



    subroutine bc_set_destructor(self)

        class(bc_set), intent(inout) :: self
        integer :: i

        do i = 1, self%m_n_bcs
            call bc_destructor_wrapper(self%m_bcs(i)%content)
            deallocate(self%m_bcs(i)%content)
        enddo

        deallocate(self%m_bcs)

    end subroutine bc_set_destructor



end module bc_set_mod
