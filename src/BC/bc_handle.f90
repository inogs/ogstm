module bc_handle_mod

    use TIME_MANAGER
    use myalloc
    use bc_mod
    use rivers_mod
    use sponge_mod
    use hard_open_mod
    use nudging_mod

    implicit none

contains



    !> Just a parser for the string contained in 'boundaries.nml'
    subroutine bc_string_parser(bc_string, bc_name, bc_type, namelist_file, filenames_list, periodic, nudged)

        character(len=47), intent(in) :: bc_string
        character(len=3), intent(out) :: bc_name
        character(len=1), intent(out) :: bc_type
        character(len=7), intent(out) :: namelist_file
        character(len=22), intent(out) :: filenames_list
        logical, intent(out) :: periodic
        logical, intent(out) :: nudged
        character(len=3) :: bc_type_str
        character(len=1) :: periodic_str
        character(len=1) :: nudged_str

        character(len=3), dimension(4), parameter :: avail_bc_types = (/ "RIV", "CLO", "SPO", "OPE" /)

        if (len(trim(bc_string)) < 47) then
            write(*, *) 'ERROR: boundary conditions namelist: invalid line'
            stop
        endif

        bc_name = bc_string(1:3)

        bc_type_str = bc_string(6:8)
        if (all(avail_bc_types /= bc_type_str)) then
            write(*, *) 'ERROR: boundary conditions namelist: invalid type'
            stop
        endif
        bc_type = bc_type_str(1:1)

        namelist_file = bc_string(11:17)
        if (namelist_file /= bc_name//".nml") then
            write(*, *) 'WARN: boundary conditions namelist: not a standard namelist: ', namelist_file
        endif

        filenames_list = bc_string(20:41)
        if (filenames_list /= "files_namelist_"//bc_name//".dat") then
            write(*, *) 'WARN: boundary conditions namelist: not a standard files namelist: ', filenames_list
        endif

        periodic_str = bc_string(44:44)
        if (periodic_str == "T") then
            periodic = .true.
        elseif (periodic_str == "F") then
            periodic = .false.
        else
            write(*, *) 'ERROR: boundary conditions namelist: invalid value for periodic attribute'
            stop
        endif

        nudged_str = bc_string(47:47)
        if (nudged_str == "T") then
            nudged = .true.
        elseif (nudged_str == "F") then
            nudged = .false.
        else
            write(*, *) 'ERROR: boundary conditions namelist: invalid value for nudging attribute'
            stop
        endif

    end subroutine bc_string_parser



    !> Boundary conditions factory method
    function bc_init(bc_string) result(bc_iter)

        character(len=47), intent(in) :: bc_string
        class(bc), pointer :: bc_iter
        class(bc), pointer :: bc_aux
        character(len=3) :: bc_name
        character(len=1) :: bc_type
        character(len=7) :: namelist_file
        character(len=22) :: filenames_list
        logical :: periodic
        logical :: nudged

        ! internal bc pointers are nullified every time the function is called
        ! ideally, one should have a single pointer of class bc,
        ! the issue with this is that function like constructors cannot be used;
        ! only standard methods can be used and a separate 'init' method should be provided.
        ! Anyway this is not a big deal, since it requires only the declaration of one pointer per class;
        ! the case construct was still necessary.
        type(rivers), pointer :: m_bc_rivers => null()
        type(sponge), pointer :: m_bc_sponge => null()
        type(hard_open), pointer :: m_bc_hard_open => null()
        type(nudging), pointer :: m_bc_nudging => null()

        bc_iter => null()
        bc_aux => null()

        call bc_string_parser(bc_string, bc_name, bc_type, namelist_file, filenames_list, periodic, nudged)

        ! select right nconstructor depending on bc type
        select case(bc_type)

            case("R") ! rivers

                allocate(m_bc_rivers)
                !$acc enter data create(m_bc_rivers)

                if (periodic) then
                    m_bc_rivers = rivers(bc_name, namelist_file, filenames_list, DATESTART, DATE__END)
                else
                    m_bc_rivers = rivers(bc_name, namelist_file, filenames_list)
                endif

                if (nudged) then
                    bc_aux => m_bc_rivers
                    allocate(m_bc_nudging)
                    m_bc_nudging = nudging(bc_aux, namelist_file, jptra) ! jptra is a global variable
                    bc_iter => m_bc_nudging
                else
                    bc_iter => m_bc_rivers
                endif

            case("C") ! closed

                write(*, *) 'ERROR: closed boundary not initialized, feature is not fully supported yet'
                stop

            case("S") ! sponge

                allocate(m_bc_sponge)

                if (periodic) then
                    m_bc_sponge = sponge(bc_name, namelist_file, filenames_list, DATESTART, DATE__END)
                else
                    m_bc_sponge = sponge(bc_name, namelist_file, filenames_list)
                endif

                if (nudged) then
                    bc_aux => m_bc_sponge
                    allocate(m_bc_nudging)
                    m_bc_nudging = nudging(bc_aux, namelist_file, jptra) ! jptra is a global variable
                    bc_iter => m_bc_nudging
                else
                    bc_iter => m_bc_sponge
                endif

            case("O") ! hard-open

                allocate(m_bc_hard_open)
                !$acc enter data create(m_bc_hard_open)

                if (periodic) then
                    m_bc_hard_open = hard_open(bc_name, namelist_file, filenames_list, jptra, DATESTART, DATE__END) ! jptra is a global variable
                else
                    m_bc_hard_open = hard_open(bc_name, namelist_file, filenames_list, jptra) ! jptra is a global variable
                endif

                if (nudged) then
                    bc_aux => m_bc_hard_open
                    allocate(m_bc_nudging)
                    m_bc_nudging = nudging(bc_aux, namelist_file, jptra) ! jptra is a global variable
                    bc_iter => m_bc_nudging
                else
                    bc_iter => m_bc_hard_open
                endif

        end select

    end function bc_init



    subroutine bc_destructor_wrapper(bc_iter)

        class(bc), target, intent(inout) :: bc_iter
        class(bc), pointer :: m_bc => null()

        m_bc => bc_iter

        select type (m_bc)
            type is (bc)
                call m_bc%bc_destructor()
            class is (rivers)
                call m_bc%rivers_destructor()
            class is (sponge)
                call m_bc%sponge_destructor()
            class is (hard_open)
                call m_bc%hard_open_destructor()
            class is (nudging)
                call m_bc%nudging_destructor()
        end select

    end subroutine bc_destructor_wrapper



    ! TO DO: double check if the original flow is correctly reproduced
    subroutine bc_update(bc_iter, current_time_string)

        class(bc), target, intent(inout) :: bc_iter
        character(len=17), intent(in) :: current_time_string
        double precision :: weight
        logical :: new_data

        class(bc), pointer :: m_bc => null()

        ! write(*, *) 'INFO: inside update_bc'

        m_bc => bc_iter

        ! write(*, *) 'INFO: m_bc associated'

        select case(nsptint)

        case(0) ! no time interpolation

            ! new data is true if current time belongs to a new interval
            call m_bc%set_current_interval(current_time_string, new_data)
            ! write(*, *) 'INFO: successfully called set_current_interval'
            
            if (current_time_string == DATESTART) then
                call m_bc%load(m_bc%get_prev_idx())
                ! write(*, *) 'INFO: successfully called load'
                call m_bc%swap()
                ! write(*, *) 'INFO: successfully called swap'
                call m_bc%load(m_bc%get_next_idx())
                ! write(*, *) 'INFO: successfully called load'
#ifdef _OPENACC
                call m_bc%update_device()
#endif
                call m_bc%actualize(1.0d0) ! weight = 1.0
                ! write(*, *) 'INFO: successfully called actualize'
            else if (new_data) then
                call m_bc%swap()
                ! write(*, *) 'INFO: successfully called swap'
                call m_bc%load(m_bc%get_next_idx())
                ! write(*, *) 'INFO: successfully called load'
#ifdef _OPENACC
                call m_bc%update_device()
#endif
                call m_bc%actualize(1.0d0) ! weight = 1.0
                ! write(*, *) 'INFO: successfully called actualize'
            endif

        case(1) ! time interpolation

            ! new data is true if interpolation is occurring on a new interval
            weight = m_bc%get_interpolation_factor(current_time_string, new_data)
            ! write(*, *) 'INFO: successfully called get_interpolation_factor'
            
            if (current_time_string == DATESTART) then
                call m_bc%load(m_bc%get_prev_idx())
                ! write(*, *) 'INFO: successfully called load'
                call m_bc%swap()
                ! write(*, *) 'INFO: successfully called swap'
                call m_bc%load(m_bc%get_next_idx())
                ! write(*, *) 'INFO: successfully called load'
#ifdef _OPENACC
                call m_bc%update_device()
#endif
            else if (new_data) then
                call m_bc%swap()
                ! write(*, *) 'INFO: successfully called swap'
                call m_bc%load(m_bc%get_next_idx())
                ! write(*, *) 'INFO: successfully called load'
#ifdef _OPENACC
                call m_bc%update_device()
#endif
            endif

            call m_bc%actualize(weight)
            ! write(*, *) 'INFO: successfully called actualize'

        end select
        
    end subroutine bc_update



end module bc_handle_mod
