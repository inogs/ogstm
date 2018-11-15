!> Auxiliary subroutines and functions
module bc_aux_mod

    implicit none

    contains



        !> This is exactly the definition of 'handle_err1' which is provided in 'IOnc.f90'.

        !> The only reason why it is copied here is that the definition is not inside a module.
        subroutine handle_err1(status, mycount, fileNetCDF)
            
            use netcdf
            
            integer status, mycount
            character fileNetCDF*(*)
            
            mycount = mycount + 1
            if(status .ne. nf90_NoErr) then
                write(*,*) 'netcdf call', mycount, 'with status = ', status
                write(*,*) 'file :', fileNetCDF
                write(*,*) nf90_strerror(status)
                write(*,*) 'Stopped'
                stop 1
            endif

        end subroutine handle_err1



        !> This is exactly the definition of 'handle_err2' which is provided in 'IOnc.f90'.

        !> The only reason why it is copied here is that the definition is not inside a module.
        subroutine handle_err2(status, fileNetCDF, varname)
            
            use netcdf
            
            integer status
            character fileNetCDF*(*), varname*(*)
            
            if(status .ne. nf90_NoErr) then
                write(*,*) 'ERROR in Var = ', varname, ' file :', fileNetCDF
            endif

        end subroutine handle_err2



        !> This is exactly the definition of 'getDIMENSION' which is provided in 'IOnc.f90'.

        !> The only reason why it is copied here is that the definition is not inside a module.
        subroutine getDIMENSION(fileNetCDF, dimname, n)
            
            use netcdf
            
            implicit none
            
            character, intent(in) :: fileNetCDF*(*), dimname*(*)
            integer, intent(inout) :: n
            
            ! local
            integer DIMid, ncid, stat
            character(len=100) junk
            integer counter
            
            counter = 0
            stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)
            call handle_err1(stat, counter,FileNetCDF)
            stat = nf90_inq_dimid (ncid, dimname, DIMid)
            call handle_err1(stat, counter,FileNetCDF)
            stat = nf90_Inquire_Dimension (ncid, DIMid, junk, n)
            call handle_err1(stat, counter,FileNetCDF)
            stat = nf90_close(ncid)
            call handle_err1(stat, counter,FileNetCDF)
        
        end subroutine getDIMENSION



        !> This is exactly the definition of 'readnc_int_1d' which is provided in 'IOnc.f90'.

        !> The only reason why it is copied here is that the definition is not inside a module.
        subroutine readnc_int_1d(fileNetCDF, varname, dim1, ARRAY)
            
            use netcdf
            ! use myalloc ! included in original version, but useless
            
            implicit none
            
            character, intent(in) :: fileNetCDF*(*), varname*(*)
            integer, intent(in) :: dim1
            integer, intent(inout), dimension(dim1) :: ARRAY
            
            integer ncid, stat, VARid
            integer counter
            
            counter = 0
            
            stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)
            call handle_err1(stat, counter, FileNetCDF)
            stat = nf90_inq_varid(ncid, varname, VARid)
            call handle_err2(stat, fileNetCDF, varname)
            call handle_err1(stat, counter, FileNetCDF)
            stat = nf90_get_var(ncid, VARid, ARRAY)
            call handle_err2(stat, fileNetCDF, varname)
            call handle_err1(stat, counter, FileNetCDF)
            stat = nf90_close(ncid)
            call handle_err1(stat, counter, FileNetCDF)

        end subroutine readnc_int_1d



        !> This is exactly the definition of 'readnc_double_1d' which is provided in 'IOnc.f90'.

        !> The only reason why it is copied here is that the definition is not inside a module.
        subroutine readnc_double_1d(fileNetCDF, varname, dim1, ARRAY)

            use netcdf
            ! use myalloc ! included in original version, but useless

            implicit none

            character, intent(in) :: fileNetCDF*(*), varname*(*)
            integer, intent(in) :: dim1
            double precision, intent(inout), dimension(dim1) :: ARRAY

            integer ncid, stat, VARid
            integer counter

            counter=0

            stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)
            call handle_err1(stat, counter, FileNetCDF)
            stat = nf90_inq_varid(ncid, varname, VARid)
            call handle_err2(stat, fileNetCDF, varname)
            call handle_err1(stat, counter, FileNetCDF)
            stat = nf90_get_var(ncid, VARid, ARRAY)
            call handle_err2(stat, fileNetCDF, varname)
            call handle_err1(stat, counter, FileNetCDF)
            stat = nf90_close(ncid)
            call handle_err1(stat, counter, FileNetCDF)

        end subroutine readnc_double_1d
        
        
        
        !> This is exactly the definition of 'readnc_slice_float_2d' which is provided in 'IOnc.f90'.

        !> The only reason why it is copied here is that the definition is not inside a module.
        !> and copy_in is in double precision
        subroutine readnc_slice_float_2d(fileNetCDF, varname, M, shift)
         
            use netcdf
            use myalloc, only: nimpp, njmpp
            use modul_param, only: jpk, jpj, jpi

            implicit none
            
            character, intent(in) :: fileNetCDF*(*), varname*(*)
            integer, intent(in) :: shift
            double precision, intent(inout) :: M(jpj,jpi)
            double precision, allocatable, dimension(:,:) :: copy_in
            
            integer ncid, stat, VARid, i, j, k
            integer counter
            integer thecount(3), start(3)
            
            allocate(copy_in(jpi,jpj))
            counter = 0
            start = (/nimpp + shift, njmpp, 1/)
            thecount = (/jpi, jpj, 1/)
            
            stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)
            call handle_err1(stat, counter, FileNetCDF)
            stat = nf90_inq_varid (ncid, varname, VARid)
            call handle_err2(stat, fileNetCDF, varname)
            call handle_err1(stat, counter, FileNetCDF)
            stat = nf90_get_var(ncid, VARid, copy_in, start, thecount)
            call handle_err2(stat, fileNetCDF, varname)
            call handle_err1(stat, counter, FileNetCDF)
            stat = nf90_close(ncid)
            call handle_err1(stat, counter, FileNetCDF)
            
            do i=1, jpi
                do j=1, jpj
                    M(j,i) = copy_in(i,j)
                enddo
            enddo
            
            deallocate(copy_in)
            
        end subroutine readnc_slice_float_2d



        !> This is exactly the definition of 'readnc_slice_float' which is provided in 'IOnc.f90'.

        !> The only reason why it is copied here is that the definition is not inside a module.
        subroutine readnc_slice_float(fileNetCDF, varname, M, shift)
            
            use netcdf
            use myalloc, only: nimpp, njmpp
            use modul_param, only: jpk, jpj, jpi
            
            implicit none
            
            character,intent(in) :: fileNetCDF*(*), varname*(*)
            integer, intent(in) :: shift
            double precision, intent(inout) :: M(jpk, jpj, jpi)
            
            real, allocatable, dimension(:, :, :) :: copy_in
            integer ncid, stat, VARid, i, j, k
            integer counter
            integer thecount(4), start(4)
            
            allocate(copy_in(jpi, jpj, jpk))
            counter = 0
            start = (/ nimpp + shift, njmpp, 1, 1 /)
            thecount = (/ jpi, jpj, jpk, 1 /)
            
            stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)
            call handle_err1(stat, counter, FileNetCDF)
            stat = nf90_inq_varid(ncid, varname, VARid)
            call handle_err2(stat, fileNetCDF, varname)
            call handle_err1(stat, counter, FileNetCDF)
            stat = nf90_get_var(ncid, VARid, copy_in, start, thecount)
            call handle_err2(stat, fileNetCDF, varname)
            call handle_err1(stat, counter, FileNetCDF)
            stat = nf90_close(ncid)
            call handle_err1(stat, counter, FileNetCDF)
            
            do i=1, jpi
                do j=1, jpj
                    do k=1, jpk
                        M(k, j, i) = real(copy_in(i, j, k), 8)
                    enddo
                enddo
            enddo
            
            deallocate(copy_in)
        
        end subroutine readnc_slice_float



        !> This is exactly the definition of 'COUNT_InSubDomain' which is provided in 'domrea.f90'.

        !> The only 2 reasons why it is copied here are:
        !! - the definition is not inside a module;
        !! - since it is requiring global variables, it is preferrable to have their dependencies outside the class.
        integer function COUNT_InSubDomain_2d(sizeGLO, idxtGLOBAL)

            use modul_param, only: jpk, jpj, jpi
            use myalloc, only: idxt, tmask ! added tmask

            implicit none

            integer, intent(in) :: sizeGLO
            integer, intent(in) :: idxtGLOBAL(sizeGLO)

            ! local
            integer kk, jj, ii, jv
            integer counter, junk

            counter = 0
            do ii = 1, jpi
                do jj = 1, jpj
                    do kk = 1, jpk
                        if (tmask(kk, jj, ii) .eq. 1) then
                            junk = idxt(kk, jj, ii)
                            do jv = 1, sizeGLO
                                if (junk .eq. idxtGLOBAL(jv)) then
                                    counter = counter + 1
                                endif
                            enddo
                        endif
                    enddo
                enddo
            enddo

            COUNT_InSubDomain_2d = counter

        end function COUNT_InSubDomain_2d



        !> This is exactly the definition of 'COUNT_InSubDomain_GIB' which is provided in 'domrea.f90'.

        !> The only 2 reasons why it is copied here are:
        !! - the definition is not inside a module;
        !! - since it is requiring global variables, it is preferrable to have their dependencies outside the class.
        integer function COUNT_InSubDomain_3d(sizeGLO, idxtGLOBAL)

            use modul_param, only: jpk, jpj, jpi
            use myalloc, only: idxt, tmask ! added tmask

            implicit none

            integer, intent(in) :: sizeGLO
            integer, intent(in) :: idxtGLOBAL(sizeGLO)

            ! local
            integer kk, jj, ii, jv
            integer counter, junk

            counter = 0
            do ii = 1, jpi
                do jj = 1, jpj
                    do kk = 1, jpk
                        if (tmask(kk, jj, ii) .eq. 1) then
                            junk = idxt(kk, jj, ii)
                            do jv = 1, sizeGLO
                                if (junk .eq. idxtGLOBAL(jv)) then
                                    counter = counter + 1
                                    exit
                                endif
                            enddo
                        endif
                    enddo
                enddo
            enddo

            COUNT_InSubDomain_3d = counter

        end function COUNT_InSubDomain_3d



        !> This is exactly the definition of 'RIVRE_Indexing' which is provided in 'domrea.f90'.

        !> The only 2 reasons why it is copied here are:
        !! - the definition is not inside a module;
        !! - since it is requiring global variables, it is preferrable to have their dependencies outside the class.
        !! Notes:
        !! - it has been necessary to pass class members as function arguments;
        !! - function has been mapped into subroutine, since return value is not used.
        subroutine RE_Indexing_2d(sizeglo, idxtglo, sizeloc, ridxt)

            use modul_param, only: jpk, jpj, jpi
            use myalloc, only: idxt

            implicit none

            integer(4), intent(in) :: sizeglo
            integer(4), intent(in) :: idxtglo(sizeglo)
            integer(4), intent(in) :: sizeloc
            integer(4), intent(out) :: ridxt(4, sizeloc)

            ! local
            integer jj, ii, jv
            integer counter, junk
            
            counter = 0
            do ii = 1, jpi
                do jj = 1, jpj
                    junk = idxt(1, jj, ii)
                    do jv = 1, sizeglo
                        if (junk .eq. idxtglo(jv)) then
                            counter = counter + 1
                            ridxt(1, counter) = jv
                            ridxt(2, counter) = 1
                            ridxt(3, counter) = jj
                            ridxt(4, counter) = ii
                        endif
                    enddo
                enddo
            enddo
            
        end subroutine RE_Indexing_2d



        !> This is exactly the definition of 'GIBRE_Indexing' which is provided in 'domrea.f90'.

        !> The only 2 reasons why it is copied here are:
        !! - the definition is not inside a module;
        !! - since it is requiring global variables, it is preferrable to have their dependencies outside the class.
        !! Notes:
        !! - it has been necessary to pass class members as function arguments;
        !! - function has been mapped into subroutine, since return value is not used.
        subroutine RE_Indexing_3d(sizeglo, idxtglo, sizeloc, ridxt)

            use modul_param, only: jpk, jpj, jpi
            use myalloc, only: idxt

            implicit none

            integer(4), intent(in) :: sizeglo
            integer(4), intent(in) :: idxtglo(sizeglo)
            integer(4), intent(in) :: sizeloc
            integer(4), intent(out) :: ridxt(4, sizeloc)

            ! local
            integer kk, jj, ii, jv
            integer counter, junk
            
            counter = 0
            do ii = 1, jpi
                do jj = 1, jpj
                    do kk = 1, jpk
                        junk = idxt(kk, jj, ii)
                        do jv = 1, sizeglo
                            if (junk .eq. idxtglo(jv)) then
                                counter = counter + 1
                                ridxt(1, counter) = jv
                                ridxt(2, counter) = kk
                                ridxt(3, counter) = jj
                                ridxt(4, counter) = ii
                            endif
                        enddo
                    enddo
                enddo
            enddo
            
        end subroutine RE_Indexing_3d



end module bc_aux_mod
