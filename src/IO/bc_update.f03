module bc_update_mod

    use TIME_MANAGER
    use myalloc
    use bc_mod

    implicit none

    class(bc), pointer :: m_bc => null()

contains

    ! TO DO: double check if the original flow is correctly reproduced
    subroutine update_bc(bc_iter, current_time_string)

        class(bc), target, intent(inout) :: bc_iter
        character(len=17), intent(in) :: current_time_string
        double precision :: weight
        logical :: new_data

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
                call m_bc%actualize(1.0d0) ! weight = 1.0
                ! write(*, *) 'INFO: successfully called actualize'
            else if (new_data) then
                call m_bc%swap()
                ! write(*, *) 'INFO: successfully called swap'
                call m_bc%load(m_bc%get_next_idx())
                ! write(*, *) 'INFO: successfully called load'
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
            else if (new_data) then
                call m_bc%swap()
                ! write(*, *) 'INFO: successfully called swap'
                call m_bc%load(m_bc%get_next_idx())
                ! write(*, *) 'INFO: successfully called load'
            endif

            call m_bc%actualize(weight)
            ! write(*, *) 'INFO: successfully called actualize'

        end select
        
    end subroutine update_bc

end module bc_update_mod
