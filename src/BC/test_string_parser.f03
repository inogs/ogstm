module functions_to_be_tested

    implicit none

contains

    subroutine bc_string_parser(bc_string, bc_name, bc_type, namelist_file, filenames_list, periodic, nudging)

        character(len=47), intent(in) :: bc_string

        character(len=3), intent(out) :: bc_name
        character(len=1), intent(out) :: bc_type
        character(len=7), intent(out) :: namelist_file
        character(len=22), intent(out) :: filenames_list
        logical, intent(out) :: periodic
        logical, intent(out) :: nudging

        character(len=3) :: bc_type_str
        character(len=1) :: periodic_str
        character(len=1) :: nudging_str

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

        nudging_str = bc_string(47:47)
        if (nudging_str == "T") then
            nudging = .true.
        elseif (nudging_str == "F") then
            nudging = .false.
        else
            write(*, *) 'ERROR: boundary conditions namelist: invalid value for nudging attribute'
            stop
        endif

    end subroutine bc_string_parser

end module functions_to_be_tested



program test_bc_string_parser

    use functions_to_be_tested

    implicit none

    character(len=47), parameter :: right_str = "riv, RIV, riv.nml, files_namelist_riv.dat, T, T"
    !character(len=47), parameter :: wrong_str_1 = "riv, RIV,  riv.nml, files_namelist_riv.dat, T, T"
    !character(len=47), parameter :: wrong_str_2 = "riv,  RIV, riv.nml, files_namelist_riv.dat, T,x"
    !character(len=47), parameter :: wrong_str_3 = "riv, RIV, riv.nml,  files_namelist_riv.dat, T,x"
    !character(len=47), parameter :: wrong_str_4 = "riv, RIV, riv.nml, files_namelist_xxx.dat, T, T"
    !character(len=47), parameter :: wrong_str_5 = "riv, RIV, riv.nml, files_namelist_riv.dat, Y, T"

    character(len=3) :: bc_name
    character(len=1) :: bc_type
    character(len=7) :: namelist_file
    character(len=22) :: filenames_list
    logical :: periodic
    logical :: nudging

    call bc_string_parser(right_str, bc_name, bc_type, namelist_file, filenames_list, periodic, nudging)
    !call bc_string_parser(wrong_str_1, bc_name, bc_type, namelist_file, filenames_list, periodic, nudging)
    !call bc_string_parser(wrong_str_2, bc_name, bc_type, namelist_file, filenames_list, periodic, nudging)
    !call bc_string_parser(wrong_str_3, bc_name, bc_type, namelist_file, filenames_list, periodic, nudging)
    !call bc_string_parser(wrong_str_4, bc_name, bc_type, namelist_file, filenames_list, periodic, nudging)
    !call bc_string_parser(wrong_str_5, bc_name, bc_type, namelist_file, filenames_list, periodic, nudging)

    write(*, *) 'bc_name = ', bc_name
    write(*, *) 'bc_type = ', bc_type
    write(*, *) 'namelist_file = ', namelist_file
    write(*, *) 'filenames_list = ', filenames_list
    write(*, *) 'periodic = ', periodic
    write(*, *) 'nudging = ', nudging

end program test_bc_string_parser
