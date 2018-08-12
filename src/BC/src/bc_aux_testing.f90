module bc_aux_mod

    implicit none

    contains
        
        ! This is exactly the definition of 'handle_err1' which is provided in 'IOnc.f90'.
        ! The only reason why it is copied here is that the definition is not inside a module.
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

        ! This is exactly the definition of 'handle_err2' which is provided in 'IOnc.f90'.
        ! The only reason why it is copied here is that the definition is not inside a module.
        subroutine handle_err2(status, fileNetCDF, varname)
            
            use netcdf
            
            integer status
            character fileNetCDF*(*), varname*(*)
            
            if(status .ne. nf90_NoErr) then
                write(*,*) 'ERROR in Var = ', varname, ' file :', fileNetCDF
            endif

        end subroutine handle_err2

        ! This is exactly the definition of 'getDIMENSION' which is provided in 'IOnc.f90'.
        ! The only reason why it is copied here is that the definition is not inside a module.
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

        ! This is exactly the definition of 'readnc_int_1d' which is provided in 'IOnc.f90'.
        ! The only reason why it is copied here is that the definition is not inside a module.
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

        ! WARNING: this is not the actual 'COUNT_InSubDomain_GIB' function,
        ! but just a replacement in order to perform serial unit testing on sponge class.
        ! The actual version of the function is the commented one above.
        ! TO DO: this should be avoided and full mpi tests enabled.
        integer(4) function COUNT_InSubDomain_GIB(sizeGLO, idxtGLOBAL)

            integer, intent(in) :: sizeGLO
            integer, intent(in) :: idxtGLOBAL(sizeGLO)

            COUNT_InSubDomain_GIB = 1

        end function COUNT_InSubDomain_GIB

end module bc_aux_mod
