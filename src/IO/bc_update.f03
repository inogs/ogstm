module bc_update_mod

    use TIME_MANAGER
    use bc_mod

    implicit none

    class(bc), pointer :: m_bc => null()

contains

    subroutine update_bc(bc_iter, current_time_string)

        class(bc), target, intent(inout) :: bc_iter
        character(len=17), intent(in) :: current_time_string
        double precision :: weight
        logical :: new_data

        m_bc => bc_iter
        ! new data is true if interpolation is occurring on a new interval
        weight = m_bc%get_interpolation_factor(current_time_string, new_data)

        ! corner cases
        if (current_time_string == DATESTART) then
            call m_bc%load(m_bc%get_prev_idx())
            call m_bc%swap()
            call m_bc%load(m_bc%get_next_idx())
        else if (new_data) ! Should it be a new if?
            call m_bc%swap()
            call m_bc%load(m_bc%get_next_idx())
        endif

        ! write 'set_indexes' method both in bc and bc_data class
        ! and apply the switch between time interpolation yes / no from the beginning

    end subroutine update_bc

end module bc_update_mod
