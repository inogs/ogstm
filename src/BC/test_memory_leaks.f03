program test_memory_leaks

    use bc_data_mod

    implicit none

    type(bc_data), pointer :: m_bc_data => null()
    character(len=24), dimension(4) :: reference_list
    character(len=24) :: tmp_file
    integer :: i ! counter

    reference_list = [ &
        "GIB_20170215-12:00:00.nc", &
        "GIB_20170515-12:00:00.nc", &
        "GIB_20170815-12:00:00.nc", &
        "GIB_20171115-12:00:00.nc" &
        ]

    allocate(m_bc_data)
    m_bc_data = bc_data("files_namelist.dat")

    do i = 1, 4
        tmp_file = m_bc_data%get_file_by_index(i)
        if (tmp_file .eq. reference_list(i)) then
            write(*, *) 'PASSED'
        else
            write(*, *) 'FAILED'
        endif
    enddo

    deallocate(m_bc_data)
    write(*, *) 'INFO: m_bc_data deallocated'
    nullify(m_bc_data)
    write(*, *) 'INFO: m_bc_data deassociated'

end program test_memory_leaks
