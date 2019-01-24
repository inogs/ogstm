program test_memory_leaks

    use bc_mod
    use sponge_mod
    use nudging_mod

    implicit none

    type(sponge), pointer :: gibraltar_no_nudging => null()
    type(nudging), pointer :: gibraltar_nudging => null()
    class(bc), pointer :: gibraltar => null()
    character(len=27), parameter :: reference = "BC/GIB_20161115-12:00:00.nc"

    allocate(gibraltar_no_nudging)
    gibraltar_no_nudging = sponge("gib", "gib.nml", "files_namelist_gib.dat")
    allocate(gibraltar_nudging)
    gibraltar_nudging = nudging(gibraltar_no_nudging, "gib.nml", 51)
    gibraltar => gibraltar_nudging

    ! do something
    if (gibraltar%get_file_by_index(1) == reference) then
        write(*, *) 'INFO: gibraltar initialized correctly'
    else
        write(*, *) 'WARN: gibraltar initialization failed'
    endif

    ! deallocate
    select type(gibraltar)
        class is (nudging)
            call gibraltar%nudging_destructor()
    end select

end program test_memory_leaks
