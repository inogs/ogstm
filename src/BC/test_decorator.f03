program test_decorator

    use nudging_mod
    use sponge_mod

    implicit none

    type(sponge), pointer :: sponge_no_nudging => null()
    type(nudging), pointer :: sponge_nudging => null()

    allocate(sponge_no_nudging)
    sponge_no_nudging = sponge( &
        "files_namelist_gib.dat", &
        "gib", &
        7, &
        "O2o N1p N3n N5s O3c O3h N6r" &
    )

    allocate(sponge_nudging)
    sponge_nudging = nudging("files_namelist_gib.dat", sponge_no_nudging, 1)

    call sponge_nudging%load(1)

    call sponge_no_nudging%sponge_destructor()
    deallocate(sponge_no_nudging)
    write(*, *) 'INFO: sponge_no_nudging deallocated'
    nullify(sponge_no_nudging)
    write(*, *) 'INFO: sponge_no_nudging deassociated'

    call sponge_nudging%nudging_destructor()
    deallocate(sponge_nudging)
    write(*, *) 'INFO: sponge_nudging deallocated'
    nullify(sponge_nudging)
    write(*, *) 'INFO: sponge_nudging deassociated'

end program test_decorator
