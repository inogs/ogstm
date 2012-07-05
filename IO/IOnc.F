       SUBROUTINE EXISTVAR(fileNetCDF, varname,B)
       USE netcdf
       IMPLICIT NONE
       character, INTENT(IN) :: fileNetCDF*(*) ,varname*(*)
       LOGICAL,   INTENT(OUT)::B
       ! local
       integer stat,ncid,VARid,counter

       counter=0
       B=.false.
         stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)  ; call handle_err1(stat, counter,FileNetCDF)
         stat = nf90_inq_varid (ncid, varname, VARid)      ;

         if(stat .ne. nf90_NoErr)  then
           !write(*,*) 'ERROR in Var = ', varname, ' file :', fileNetCDF
           B = .false.
         else
            B=.true.
        endif

       END SUBROUTINE EXISTVAR

! *********************************************************************************

      SUBROUTINE readnc_slice_double(fileNetCDF,varname, M)
      USE myalloc
      USE myalloc_mpp
      USE netcdf
      implicit none


      character fileNetCDF*(*) ,varname*(*)
      integer ncid, stat, VARid

      integer counter
      integer thecount(4), start(4)
      real(8) M(jpi,jpj,jpk)

      counter = 0;
      start    = (/nimpp, njmpp,  1,  1/)
      thecount = (/jpi,     jpj, jpk, 1/)

      stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)  ; call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_inq_varid (ncid, varname, VARid)      ; call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_get_var (ncid,VARid,M,start, thecount);
      call handle_err2(stat, fileNetCDF,varname)        ;  call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_close(ncid)                           ; call handle_err1(stat, counter,FileNetCDF)

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

      integer counter
      integer thecount(4), start(4)
      real(8) M(jpi,jpj,jpk)
      real(4) F(jpi,jpj,jpk)

      counter = 0;
      start    = (/nimpp, njmpp,  1,  1/)
      thecount = (/jpi,     jpj, jpk, 1/)

      stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)  ; call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_inq_varid (ncid, varname, VARid)      ; call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_get_var (ncid,VARid,F,start, thecount);
      call handle_err2(stat, fileNetCDF,varname)        ; call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_close(ncid)                           ; call handle_err1(stat, counter,FileNetCDF)

      M = real(F,8);


      END SUBROUTINE readnc_slice_float

      SUBROUTINE readnc_slice_int(fileNetCDF,varname, M)
      USE myalloc
      USE netcdf
      implicit none


      character fileNetCDF*(*) ,varname*(*)
      integer ncid, stat, VARid

      integer counter
      integer thecount(4), start(4)
      integer M(jpi,jpj,jpk)
      counter = 0;
      start    = (/nimpp, njmpp,  1,  1/)
      thecount = (/jpi,     jpj, jpk, 1/)

      stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)  ; call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_inq_varid (ncid, varname, VARid)      ; call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_get_var (ncid,VARid,M,start, thecount)
      call handle_err2(stat, fileNetCDF,varname)        ; call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_close(ncid)                           ; call handle_err1(stat, counter,FileNetCDF)

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

      integer counter
      integer thecount(3), start(3)
      real(8) M(jpi,jpj)


      counter = 0;
      start    = (/nimpp, njmpp,  1/)
      thecount = (/jpi,     jpj,  1/)


      stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)  ; call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_inq_varid (ncid, varname, VARid)      ; call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_get_var (ncid,VARid,M,start, thecount);
      call handle_err2(stat, fileNetCDF,varname)        ; call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_close(ncid)                           ; call handle_err1(stat, counter,FileNetCDF)




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

      integer counter
      integer thecount(3), start(3)
      real(8) M(jpi,jpj)
      real(4) F(jpi,jpj)

      counter = 0;
      start    = (/nimpp, njmpp,  1/)
      thecount = (/jpi,     jpj,  1/)


      stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)  ; call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_inq_varid (ncid, varname, VARid)      ; call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_get_var (ncid,VARid,F,start, thecount)
      call handle_err2(stat, fileNetCDF,varname)        ; call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_close(ncid)                           ; call handle_err1(stat, counter,FileNetCDF)

      M = real(F,8);


      END SUBROUTINE readnc_slice_float_2d


! ************************************************************
      SUBROUTINE readnc_global_double(fileNetCDF,varname, M)
      USE myalloc
      USE netcdf
      implicit none


      character fileNetCDF*(*) ,varname*(*)
      integer ncid, stat, VARid

      integer counter
      integer thecount(4), start(4)
      real(8) M(jpi,jpj,jpk)

      counter = 0;
      start    = (/1,       1,       1,  1/)
      thecount = (/jpiglo,  jpjglo, jpk, 1/)


      stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)  ; call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_inq_varid (ncid, varname, VARid)      ; call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_get_var (ncid,VARid,M,start, thecount)
      call handle_err2(stat, fileNetCDF,varname)        ; call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_close(ncid)                           ; call handle_err1(stat, counter,FileNetCDF)

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

        stat = nf90_open(fileNetCDF, nf90_nowrite, ncid) ; call handle_err1(stat, counter,FileNetCDF)
        stat = nf90_inq_varid (ncid, varname, VARid)     ; call handle_err1(stat, counter,FileNetCDF)
        stat = nf90_get_var (ncid,VARid,ARRAY,start, thecount);
        call handle_err2(stat, fileNetCDF,varname)       ; call handle_err1(stat, counter,FileNetCDF)
        stat = nf90_close(ncid)                          ; call handle_err1(stat, counter,FileNetCDF)


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

        stat = nf90_open(fileNetCDF, nf90_nowrite, ncid) ; call handle_err1(stat, counter,FileNetCDF)
        stat = nf90_inq_varid (ncid, varname, VARid)     ; call handle_err1(stat, counter,FileNetCDF)
        stat = nf90_get_var (ncid,VARid,ARRAY)           ;
        call handle_err2(stat, fileNetCDF,varname)       ; call handle_err1(stat, counter,FileNetCDF)
        stat = nf90_close(ncid)                          ; call handle_err1(stat, counter,FileNetCDF)


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

        stat = nf90_open(fileNetCDF, nf90_nowrite, ncid) ; call handle_err1(stat, counter,FileNetCDF)
        stat = nf90_inq_varid (ncid, varname, VARid)     ; call handle_err1(stat, counter,FileNetCDF)
        stat = nf90_get_var (ncid,VARid,ARRAY)
        call handle_err2(stat, fileNetCDF,varname)       ; call handle_err1(stat, counter,FileNetCDF)
        stat = nf90_close(ncid)                          ; call handle_err1(stat, counter,FileNetCDF)


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
        stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)    ; call handle_err1(stat, counter,FileNetCDF)
        stat = nf90_inq_dimid (ncid, dimname, DIMid)        ; call handle_err1(stat, counter,FileNetCDF)
        stat = nf90_Inquire_Dimension (ncid, DIMid, junk, n); call handle_err1(stat, counter,FileNetCDF)
        stat = nf90_close(ncid)                             ; call handle_err1(stat, counter,FileNetCDF)
        END SUBROUTINE getDIMENSION

       !****************************************************************************
       SUBROUTINE write_restart(fileNetCDF,varB,varN)
       USE netcdf
       USE myalloc

       IMPLICIT NONE
       CHARACTER*(*) fileNetCDF
       real(8) varB, varN
       CHARACTER(LEN=3) VAR
       CHARACTER(LEN=6) RSTVAR
       CHARACTER(LEN=17) TimeString

      ! local
       integer s, nc, counter
       integer timid, depid, yid, xid, xaid, yaid, zaid
       integer idB, idN, idLon, idLat, idLev, idTim

       TimeString =fileNetCDF( 5:21)
       VAR        =fileNetCDF(23:25)


      s = nf90_create(fileNetCDF, NF90_CLOBBER, nc)

      s = nf90_put_att(nc, nf90_global, 'TimeString'     , TimeString)
        ! *********** DIMENSIONS ****************
        s= nf90_def_dim(nc,'x'   , jpiglo,  xid)
        s= nf90_def_dim(nc,'y'   , jpjglo,  yid)
        s= nf90_def_dim(nc,'z'   , jpk   ,depid)
        s= nf90_def_dim(nc,'time', 1     ,timid)
        s= nf90_def_dim(nc,'x_a'      , 1,  xaid)
        s= nf90_def_dim(nc,'y_a'      , 1,  yaid)
        s= nf90_def_dim(nc,'z_a'      , 3  ,zaid)


       s = nf90_def_var(nc,'nav_lon', nf90_double,  (/xid,yid/), idLon)
       s = nf90_def_var(nc,'nav_lat', nf90_double,  (/xid,yid/), idLat)
       s = nf90_def_var(nc,'nav_lev', nf90_double,  (/depid/)  , idLev)
       s = nf90_def_var(nc,'time'   , nf90_double,  (/timid/)  , idTim)


        RSTVAR='TRB'//VAR;
        s = nf90_def_var(nc,RSTVAR , nf90_double, (/xid,yid,depid,timid/),   idB)

        RSTVAR='TRN'//VAR;
        s = nf90_def_var(nc,RSTVAR , nf90_double, (/xid,yid,depid,timid/),   idN)


        s =nf90_enddef(nc)

        s =nf90_close(nc)

!      to be continued


       END SUBROUTINE write_restart

       !****************************************************************************
       ! physical data were printed as grid3D.dat grid2D.dat files
       SUBROUTINE physDump(fileNetCDF, datefrom, dateTo)
       USE myalloc
       USE netcdf

       IMPLICIT NONE

        character fileNetCDF*(*)
        character(LEN=17) datefrom, dateTo
        integer s, nc, counter
        integer timid, depid, yid, xid
        integer idvartime, idgdept, idphit, idlamt
        integer idT, idS, idU, idV, idW, idEddy, idR, idWs, idE


        s = nf90_create(fileNetCDF, NF90_CLOBBER, nc)
        ! *********** GLOBAL ********************
        s = nf90_put_att(nc, nf90_global, 'Convenctions'     , 'OPA')
        s = nf90_put_att(nc, nf90_global, 'DateStart'     , datefrom)
        s = nf90_put_att(nc, nf90_global, 'Date__End'     ,   dateTo)

        ! *********** DIMENSIONS ****************
        s= nf90_def_dim(nc,'x'           , jpiglo,  xid)
        s= nf90_def_dim(nc,'y'           , jpjglo,  yid)
        s= nf90_def_dim(nc,'deptht'      , jpk   ,depid)
        s= nf90_def_dim(nc,'time_counter', NF90_UNLIMITED,timid)


        ! ********** VARIABLES *****************
        s = nf90_def_var(nc,'time_counter',  nf90_double,(/timid/),           idvartime)
        s = nf90_def_var(nc,'deptht',        nf90_float, (/depid/),             idgdept)
        s = nf90_def_var(nc,'nav_lat',       nf90_float, (/xid,yid/),            idphit)
        s = nf90_def_var(nc,'nav_lon',       nf90_float, (/xid,yid/),            idlamt)

        s = nf90_def_var(nc,'vosaline' ,     nf90_float, (/xid,yid,depid,timid/),   idS)
        s = nf90_def_var(nc,'votemper' ,     nf90_float, (/xid,yid,depid,timid/),   idT)
        s = nf90_def_var(nc,'vozocrtx' ,     nf90_float, (/xid,yid,depid,timid/),   idU)
        s = nf90_def_var(nc,'vomecrty' ,     nf90_float, (/xid,yid,depid,timid/),   idV)
        s = nf90_def_var(nc,'vovecrtz' ,     nf90_float, (/xid,yid,depid,timid/),   idW)
        s = nf90_def_var(nc,'votkeavt' ,     nf90_float, (/xid,yid,depid,timid/),idEddy)

        s = nf90_def_var(nc,'soshfldo' ,     nf90_float, (/xid,yid,      timid/),   idR)
        s = nf90_def_var(nc,'sowindsp' ,     nf90_float, (/xid,yid,      timid/),  idWs)
        s = nf90_def_var(nc,'sowaflcd' ,     nf90_float, (/xid,yid,      timid/),   idE)


        !s = nf90_put_att(nc,idvartime,'units', UnitsTime )

        s = nf90_put_att(nc,idgdept,'units'        ,'m')
        s = nf90_put_att(nc,idgdept,'positive'     ,'down')
        s = nf90_put_att(nc,idphit, 'units'        ,'degrees_north')
        s = nf90_put_att(nc,idphit, 'long_name'    ,'Latitude')
        s = nf90_put_att(nc,idlamt, 'units'        ,'degrees_east')
        s = nf90_put_att(nc,idlamt, 'long_name'    ,'Longitude')

        s = nf90_put_att(nc,idS   , 'long_name'    ,'Salinity')
        s = nf90_put_att(nc,idS   , 'units'        ,'PSU')
        s = nf90_put_att(nc,idS   , 'missing_value',1.e+20)

        s = nf90_put_att(nc,idT   , 'long_name'    ,'Temperature')
        s = nf90_put_att(nc,idT   , 'units'        ,'C')
        s = nf90_put_att(nc,idT   , 'missing_value',1.e+20)

        s = nf90_put_att(nc,idU   , 'long_name'    ,'Zonal Current')
        s = nf90_put_att(nc,idU   , 'units'        ,'m/s')
        s = nf90_put_att(nc,idU   , 'missing_value',1.e+20)

        s = nf90_put_att(nc,idV   , 'long_name'    ,'Meridional Current')
        s = nf90_put_att(nc,idV   , 'units'        ,'m/s')
        s = nf90_put_att(nc,idV   , 'missing_value',1.e+20)

        s = nf90_put_att(nc,idW   , 'long_name'    ,'Vertical Velocity')
        s = nf90_put_att(nc,idW   , 'units'        ,'m/s')
        s = nf90_put_att(nc,idW   , 'missing_value',1.e+20)

        s = nf90_put_att(nc,idEddy, 'long_name'    ,'Vertical Eddy Diffusivity')
        s = nf90_put_att(nc,idEddy, 'units'        ,'m/s')
        s = nf90_put_att(nc,idEddy, 'missing_value',1.e+20)

!       2D
        s = nf90_put_att(nc,idR   , 'long_name'    ,'Short_Wave_Radiation')
        s = nf90_put_att(nc,idR   , 'units'        ,'W/m2')
        s = nf90_put_att(nc,idR   , 'missing_value',1.e+20)

        s = nf90_put_att(nc,idWs  , 'long_name'    ,'wind_speed')
        s = nf90_put_att(nc,idWs  , 'units'        ,'m/s')
        s = nf90_put_att(nc,idWs  , 'missing_value',1.e+20)

        s = nf90_put_att(nc,idE   , 'long_name'    ,'Concentration/Dilution Water Flux')
        s = nf90_put_att(nc,idE   , 'units'        ,'kg/m*2/s')
        s = nf90_put_att(nc,idE   , 'missing_value',1.e+20)

        s =nf90_enddef(nc)

        counter=0

        s = nf90_put_var(nc, idlamt,   REAL(glamt,4)  ); call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idphit,   REAL(gphit,4)  ); call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idgdept,  REAL(gdept,4)  ); call handle_err1(s,counter,fileNetCDF)
!      3D vars
        s = nf90_put_var(nc, idS,      REAL(totsnIO,4)  ); call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idT,      REAL(tottnIO,4)  ); call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idU,      REAL(totunIO,4)  ); call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idV,      REAL(totvnIO,4)  ); call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idW,      REAL(totwnIO,4)  ); call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idEddy,   REAL(totavtIO,4) ); call handle_err1(s,counter,fileNetCDF)

!       2D vars
        s = nf90_put_var(nc, idWs,     REAL(totvatmIO,4) ); call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idE,      REAL(totempIO,4)  ); call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idR,      REAL(totqsrIO,4)  ); call handle_err1(s,counter,fileNetCDF)



        s= nf90_close(nc)


       END SUBROUTINE physDump

       !****************************************************************************


      !****************************************************************************
        subroutine handle_err1(status,mycount, fileNetCDF)
        USE netcdf
        integer status,mycount
        character fileNetCDF*(*)
        mycount =mycount+1
        if(status .ne. nf90_NoErr)  then
           write(*,*) 'netcdf call',mycount,'with status = ',status
           write(*,*)  'file :', fileNetCDF
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

       ! ***************************************************************************



