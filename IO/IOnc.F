      SUBROUTINE readnc_slice_double(fileNetCDF,varname, M)
      USE myalloc
      USE myalloc_mpp
      USE netcdf
      implicit none


      character fileNetCDF*(*) ,varname*(*)
      integer ncid, stat, VARid

      integer mycount
      integer thecount(4), start(4)
      real(8) M(jpi,jpj,jpk)

      mycount = 0;
      start    = (/nimpp, njmpp,  1,  1/)
      thecount = (/jpi,     jpj, jpk, 1/)

      stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)  ; call handle_err1(stat, mycount)
      stat = nf90_inq_varid (ncid, varname, VARid)      ; call handle_err1(stat, mycount)
      stat = nf90_get_var (ncid,VARid,M,start, thecount);
      call handle_err2(stat, fileNetCDF,varname);  call handle_err1(stat, mycount)
      stat = nf90_close(ncid)                           ; call handle_err1(stat, mycount)

      END SUBROUTINE readnc_slice_double

! ********************************************************************************************
! ********************************************************************************************
! ********************************************************************************************

      SUBROUTINE readnc_slice_float(fileNetCDF,varname, M)
      USE myalloc
      USE netcdf
      implicit none


      character fileNetCDF*(*) ,varname*(*)
      integer ncid, stat, VARid

      integer mycount
      integer thecount(4), start(4)
      real(8) M(jpi,jpj,jpk)
      real(4) F(jpi,jpj,jpk)

      mycount = 0;
      start    = (/nimpp, njmpp,  1,  1/)
      thecount = (/jpi,     jpj, jpk, 1/)

      stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)  ; call handle_err1(stat, mycount)
      stat = nf90_inq_varid (ncid, varname, VARid)      ; call handle_err1(stat, mycount)
      stat = nf90_get_var (ncid,VARid,F,start, thecount);
      call handle_err2(stat, fileNetCDF,varname)        ; call handle_err1(stat, mycount)
      stat = nf90_close(ncid)                           ; call handle_err1(stat, mycount)

      M = real(F,8);


      END SUBROUTINE readnc_slice_float

      SUBROUTINE readnc_slice_int(fileNetCDF,varname, M)
      USE myalloc
      USE netcdf
      implicit none


      character fileNetCDF*(*) ,varname*(*)
      integer ncid, stat, VARid

      integer mycount
      integer thecount(4), start(4)
      integer M(jpi,jpj,jpk)
      mycount = 0;
      start    = (/nimpp, njmpp,  1,  1/)
      thecount = (/jpi,     jpj, jpk, 1/)

      stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)  ; call handle_err1(stat, mycount)
      stat = nf90_inq_varid (ncid, varname, VARid)      ; call handle_err1(stat, mycount)
      stat = nf90_get_var (ncid,VARid,M,start, thecount)
      call handle_err2(stat, fileNetCDF,varname)        ; call handle_err1(stat, mycount)
      stat = nf90_close(ncid)                           ; call handle_err1(stat, mycount)

      END SUBROUTINE readnc_slice_int

! ********************************************************************************************
! ********************************************************************************************
! ********************************************************************************************

      SUBROUTINE readnc_slice_double_2d(fileNetCDF,varname, M)
      USE myalloc
      USE myalloc_mpp
      USE netcdf
      implicit none


      character fileNetCDF*(*) ,varname*(*)
      integer ncid, stat, VARid

      integer mycount
      integer thecount(3), start(3)
      real(8) M(jpi,jpj)


      mycount = 0;
      start    = (/nimpp, njmpp,  1/)
      thecount = (/jpi,     jpj,  1/)


      stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)  ; call handle_err1(stat, mycount)
      stat = nf90_inq_varid (ncid, varname, VARid)      ; call handle_err1(stat, mycount)
      stat = nf90_get_var (ncid,VARid,M,start, thecount);
      call handle_err2(stat, fileNetCDF,varname)        ; call handle_err1(stat, mycount)
      stat = nf90_close(ncid)                           ; call handle_err1(stat, mycount)




      END SUBROUTINE readnc_slice_double_2d


! ********************************************************************************************
! ********************************************************************************************
! ********************************************************************************************

      SUBROUTINE readnc_slice_float_2d(fileNetCDF,varname, M)
      USE myalloc
      USE netcdf
      implicit none


      character fileNetCDF*(*) ,varname*(*)
      integer ncid, stat, VARid

      integer mycount
      integer thecount(3), start(3)
      real(8) M(jpi,jpj)
      real(4) F(jpi,jpj)

      mycount = 0;
      start    = (/nimpp, njmpp,  1/)
      thecount = (/jpi,     jpj,  1/)


      stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)  ; call handle_err1(stat, mycount)
      stat = nf90_inq_varid (ncid, varname, VARid)      ; call handle_err1(stat, mycount)
      stat = nf90_get_var (ncid,VARid,F,start, thecount)
      call handle_err2(stat, fileNetCDF,varname)        ; call handle_err1(stat, mycount)
      stat = nf90_close(ncid)                           ; call handle_err1(stat, mycount)

      M = real(F,8);


      END SUBROUTINE readnc_slice_float_2d


! ************************************************************
      SUBROUTINE readnc_global_double(fileNetCDF,varname, M)
      USE myalloc
      USE netcdf
      implicit none


      character fileNetCDF*(*) ,varname*(*)
      integer ncid, stat, VARid

      integer mycount
      integer thecount(4), start(4)
      real(8) M(jpi,jpj,jpk)

      mycount = 0;
      start    = (/1,       1,       1,  1/)
      thecount = (/jpiglo,  jpjglo, jpk, 1/)


      stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)  ; call handle_err1(stat, mycount)
      stat = nf90_inq_varid (ncid, varname, VARid)      ; call handle_err1(stat, mycount)
      stat = nf90_get_var (ncid,VARid,M,start, thecount)
      call handle_err2(stat, fileNetCDF,varname)        ; call handle_err1(stat, mycount)
      stat = nf90_close(ncid)                           ; call handle_err1(stat, mycount)

      END SUBROUTINE readnc_global_double



        !****************************************************************************
        SUBROUTINE readmask_double_1d(fileNetCDF,varname,ARRAY)
        use netcdf
        USE modul_param
        implicit none
        character fileNetCDF*(*) ,varname*(*)
        integer ncid, stat, VARid
        integer counter
        real(8) ARRAY(jpk)
        integer thecount(4), start(4)

        counter=0
        start    = (/1,       1,       1,  1/)
        thecount = (/1,       1,      jpk, 1/)

        stat = nf90_open(fileNetCDF, nf90_nowrite, ncid) ; call handle_err1(stat,counter)
        stat = nf90_inq_varid (ncid, varname, VARid)     ; call handle_err1(stat,counter)
        stat = nf90_get_var (ncid,VARid,ARRAY,start, thecount);
        call handle_err2(stat, fileNetCDF,varname)       ; call handle_err1(stat,counter)
        stat = nf90_close(ncid)                          ; call handle_err1(stat,counter)


        end SUBROUTINE readmask_double_1d


        !****************************************************************************
        SUBROUTINE readnc_double_1d(fileNetCDF,varname,im,ARRAY)
        use netcdf
        implicit none
        character fileNetCDF*(*) ,varname*(*)
        integer ncid, stat, VARid
        integer im,counter
        real(8) ARRAY(im)

        counter=0

        stat = nf90_open(fileNetCDF, nf90_nowrite, ncid) ; call handle_err1(stat,counter)
        stat = nf90_inq_varid (ncid, varname, VARid)     ; call handle_err1(stat,counter)
        stat = nf90_get_var (ncid,VARid,ARRAY)           ;
        call handle_err2(stat, fileNetCDF,varname)       ; call handle_err1(stat,counter)
        stat = nf90_close(ncid)                          ; call handle_err1(stat,counter)


        end SUBROUTINE readnc_double_1d

        !****************************************************************************
        SUBROUTINE readnc_int_1d(fileNetCDF,varname,im,ARRAY)
        use netcdf
        implicit none
        character fileNetCDF*(*) ,varname*(*)
        integer ncid, stat, VARid
        integer im, counter
        integer ARRAY(im)

        counter=0

        stat = nf90_open(fileNetCDF, nf90_nowrite, ncid) ; call handle_err1(stat,counter)
        stat = nf90_inq_varid (ncid, varname, VARid)     ; call handle_err1(stat,counter)
        stat = nf90_get_var (ncid,VARid,ARRAY)
        call handle_err2(stat, fileNetCDF,varname)       ; call handle_err1(stat,counter)
        stat = nf90_close(ncid)                          ; call handle_err1(stat,counter)


        end SUBROUTINE readnc_int_1d

       !****************************************************************************


        SUBROUTINE getDIMENSION(fileNetCDF,dimname,n)
        use netcdf
        implicit none
        character  fileNetCDF*(*), dimname*(*)
        integer n

        ! local

        integer DIMid,ncid,stat
        character(LEN=100) junk
        integer counter

        counter = 0
        stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)    ; call handle_err1(stat,counter)
        stat = nf90_inq_dimid (ncid, dimname, DIMid)     ; call handle_err1(stat,counter)
        stat = nf90_Inquire_Dimension (ncid, DIMid, junk, n); call handle_err1(stat,counter)
        stat = nf90_close(ncid)                             ; call handle_err1(stat,counter)
        END SUBROUTINE getDIMENSION

        !****************************************************************************
       SUBROUTINE write_restart(filename,varB,varN)
       USE netcdf
       USE myalloc

       IMPLICIT NONE
       CHARACTER*(*) filename
       real(8) varB, varN

!      to be continued


       END SUBROUTINE write_restart


      !****************************************************************************
        subroutine handle_err1(status,mycount)
        USE netcdf
        integer status,mycount
        mycount =mycount+1
        if(status .ne. nf90_NoErr)  then
           write(*,*) 'netcdf call',mycount,'with status = ',status
           write(*,*) nf90_strerror(status)
           write(*,*) 'Stopped'
           STOP 1
        endif

        end subroutine handle_err1


      ! **************************************************************************
        subroutine handle_err2(status,fileNetCDF,varname)
        USE netcdf
        integer status
        character fileNetCDF*(*) ,varname*(*)

        if(status .ne. nf90_NoErr)  then
           write(*,*) 'ERROR in Var = ', varname, ' file :', fileNetCDF
        endif

        end subroutine handle_err2
