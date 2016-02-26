       SUBROUTINE EXISTVAR_bfm(fileNetCDF, varname,B)
       USE netcdf
       IMPLICIT NONE
       character, INTENT(IN) :: fileNetCDF*(*) ,varname*(*)
       LOGICAL,   INTENT(OUT)::B
       ! local
       integer stat,ncid,VARid,counter

       counter=0
       B=.false.
         stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)  ; call handle_err1(stat, counter,FileNetCDF)
         stat = nf90_inq_varid (ncid, varname, VARid)

         if(stat .ne. nf90_NoErr)  then
           !write(*,*) 'ERROR in Var = ', varname, ' file :', fileNetCDF
           B = .false.
         else
            B=.true.
        endif
        stat = nf90_close(ncid)                           ; call handle_err1(stat, counter,FileNetCDF)
       END SUBROUTINE EXISTVAR_bfm


! *****************************************************************************

      SUBROUTINE readnc_slice_double_ben(fileNetCDF,varname, M)
      USE myalloc
      USE myalloc_mpp
      USE BIOBEN_mem
      USE netcdf
      implicit none


      character fileNetCDF*(*) ,varname*(*)
      integer ncid, stat, VARid

      integer counter
      integer thecount(4), start(4)
      real(8) M(jpi,jpj,jpk_ben)

      counter = 0;
      start    = (/nimpp, njmpp,  1,  1/)
      thecount = (/jpi,     jpj, jpk_ben, 1/)

      stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)  ; call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_inq_varid (ncid, varname, VARid)      ; call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_get_var (ncid,VARid,M,start, thecount);
      call handle_err2(stat, fileNetCDF,varname)        ;  call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_close(ncid)                           ; call handle_err1(stat, counter,FileNetCDF)

      END SUBROUTINE readnc_slice_double_ben

! ********************************************************************************************
! ********************************************************************************************
! ********************************************************************************************

      SUBROUTINE readnc_slice_float_ben(fileNetCDF,varname, M)
      USE myalloc
      USE BIOBEN_mem
      USE netcdf
      implicit none


      character fileNetCDF*(*) ,varname*(*)
      integer ncid, stat, VARid

      integer counter
      integer thecount(4), start(4)
      real(8) M(jpi,jpj,jpk_ben)
      real(4) F(jpi,jpj,jpk_ben)

      counter = 0;
      start    = (/nimpp, njmpp,  1,  1/)
      thecount = (/jpi,     jpj, jpk_ben, 1/)

      stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)  ; call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_inq_varid (ncid, varname, VARid)      ; call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_get_var (ncid,VARid,F,start, thecount);
      call handle_err2(stat, fileNetCDF,varname)        ; call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_close(ncid)                           ; call handle_err1(stat, counter,FileNetCDF)

      M = real(F,8);


      END SUBROUTINE readnc_slice_float_ben

      SUBROUTINE readnc_slice_int_ben(fileNetCDF,varname, M)
      USE myalloc
      USE BIOBEN_mem
      USE netcdf
      implicit none


      character fileNetCDF*(*) ,varname*(*)
      integer ncid, stat, VARid

      integer counter
      integer thecount(4), start(4)
      integer M(jpi,jpj,jpk_ben)
      counter = 0;
      start    = (/nimpp, njmpp,  1,  1/)
      thecount = (/jpi,     jpj, jpk_ben, 1/)

      stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)  ; call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_inq_varid (ncid, varname, VARid)      ; call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_get_var (ncid,VARid,M,start, thecount)
      call handle_err2(stat, fileNetCDF,varname)        ; call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_close(ncid)                           ; call handle_err1(stat, counter,FileNetCDF)

      END SUBROUTINE readnc_slice_int_ben

! ********************************************************************************************
! ********************************************************************************************
! ********************************************************************************************

      SUBROUTINE readnc_slice_double_2d_ben(fileNetCDF,varname, M)
      USE myalloc
      USE BIOBEN_mem
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




      END SUBROUTINE readnc_slice_double_2d_ben


! ********************************************************************************************
! ********************************************************************************************
! ********************************************************************************************

      SUBROUTINE readnc_slice_float_2d_ben(fileNetCDF,varname, M)
      USE myalloc
      USE BIOBEN_mem
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


      END SUBROUTINE readnc_slice_float_2d_ben


! ************************************************************
      SUBROUTINE readnc_global_double_ben(fileNetCDF,varname, M)
      USE myalloc
      USE BIOBEN_mem
      USE netcdf
      implicit none


      character fileNetCDF*(*) ,varname*(*)
      integer ncid, stat, VARid

      integer counter
      integer thecount(4), start(4)
      real(8) M(jpiglo,jpjglo,jpk_ben)

      counter = 0;
      start    = (/1,       1,       1,  1/)
      thecount = (/jpiglo,  jpjglo, jpk_ben, 1/)


      stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)  ; call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_inq_varid (ncid, varname, VARid)      ; call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_get_var (ncid,VARid,M,start, thecount)
      call handle_err2(stat, fileNetCDF,varname)        ; call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_close(ncid)                           ; call handle_err1(stat, counter,FileNetCDF)

      END SUBROUTINE readnc_global_double_ben


! ************************************************************
      SUBROUTINE readnc_global_double_2d_ben(fileNetCDF,varname, M)
      USE myalloc
      USE BIOBEN_mem
      USE netcdf
      implicit none


      character fileNetCDF*(*) ,varname*(*)
      integer ncid, stat, VARid

      integer counter
      integer thecount(4), start(4)
      real(8) M(jpiglo,jpjglo)

      counter = 0;
      start    = (/1,       1,      1, 1/)
      thecount = (/jpiglo,  jpjglo, 1, 1/)


      stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)  ; call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_inq_varid (ncid, varname, VARid)      ; call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_get_var (ncid,VARid,M,start, thecount)
      call handle_err2(stat, fileNetCDF,varname)        ; call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_close(ncid)                           ; call handle_err1(stat, counter,FileNetCDF)

      END SUBROUTINE readnc_global_double_2d_ben



        !****************************************************************************
        SUBROUTINE readmask_double_1d_ben(fileNetCDF,varname,ARRAY)
        use netcdf
        USE BIOBEN_mem
        USE modul_param
        implicit none
        character fileNetCDF*(*) ,varname*(*)
        integer ncid, stat, VARid
        integer counter
        real(8) ARRAY(jpk_ben)
        integer thecount(4), start(4)

        counter=0
        start    = (/1,       1,       1,  1/)
        thecount = (/1,       1,      jpk_ben, 1/)

        stat = nf90_open(fileNetCDF, nf90_nowrite, ncid) ; call handle_err1(stat, counter,FileNetCDF)
        stat = nf90_inq_varid (ncid, varname, VARid)     ; call handle_err1(stat, counter,FileNetCDF)
        stat = nf90_get_var (ncid,VARid,ARRAY,start, thecount);
        call handle_err2(stat, fileNetCDF,varname)       ; call handle_err1(stat, counter,FileNetCDF)
        stat = nf90_close(ncid)                          ; call handle_err1(stat, counter,FileNetCDF)


        end SUBROUTINE readmask_double_1d_ben


        !****************************************************************************
        SUBROUTINE readnc_double_1d_ben(fileNetCDF,varname,im,ARRAY)
        USE BIOBEN_mem
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


        end SUBROUTINE readnc_double_1d_ben

        !****************************************************************************
        SUBROUTINE readnc_int_1d_ben(fileNetCDF,varname,im,ARRAY)
        USE BIOBEN_mem
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


        end SUBROUTINE readnc_int_1d_ben

       !****************************************************************************


        SUBROUTINE getDIMENSION_ben(fileNetCDF,dimname,n)
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
        END SUBROUTINE getDIMENSION_ben

       !****************************************************************************
       !****************************************************************************
       !****************************************************************************

       SUBROUTINE write_restart_ben(fileNetCDF, julian)
       USE netcdf
       USE BIOBEN_mem
       USE myalloc

       IMPLICIT NONE
       CHARACTER*(*) fileNetCDF
       REAL(8) julian

       ! local
       CHARACTER(LEN=20) VAR
       CHARACTER(LEN=17) TimeString
       integer istart, iend
       integer s, nc, counter
       integer timid, depid, yid, xid, xaid, yaid, zaid
       integer idB, idN, idLon, idLat, idLev, idTim

       TimeString =fileNetCDF(14:30)
       istart=index(fileNetCDF,'00.')+3
       iend  =index(fileNetCDF,'.nc')-1
       VAR        =fileNetCDF(istart:iend)


      s = nf90_create(fileNetCDF, NF90_CLOBBER, nc)

      s = nf90_put_att(nc, nf90_global, 'TimeString'     , TimeString)
        ! *********** DIMENSIONS ****************
        s= nf90_def_dim(nc,'x'   , jpiglo,  xid)
        s= nf90_def_dim(nc,'y'   , jpjglo,  yid)
        s= nf90_def_dim(nc,'z'   , jpk_ben   ,depid)
        s= nf90_def_dim(nc,'time', 1     ,timid)
        s= nf90_def_dim(nc,'x_a'      , 1,  xaid)
        s= nf90_def_dim(nc,'y_a'      , 1,  yaid)
        s= nf90_def_dim(nc,'z_a'      , 3  ,zaid)


       s = nf90_def_var(nc,'nav_lon', nf90_double,  (/xid,yid/), idLon)
       s = nf90_def_var(nc,'nav_lat', nf90_double,  (/xid,yid/), idLat)
       s = nf90_def_var(nc,'nav_lev', nf90_double,  (/depid/)  , idLev)
       s = nf90_def_var(nc,'time'   , nf90_double,  (/timid/)  , idTim)

        s = nf90_def_var(nc,'TRB'//trim(VAR), nf90_double, (/xid,yid,depid,timid/), idB)
        s = nf90_def_var(nc,'TRN'//trim(VAR), nf90_double, (/xid,yid,depid,timid/), idN)

        s= nf90_put_att(nc,idTim ,'Units', 'seconds since 1582-10-15 00:00:00');
        s = nf90_put_att(nc,idB   , 'missing_value',1.e+20)
        s = nf90_put_att(nc,idN   , 'missing_value',1.e+20)
        s =nf90_enddef(nc)

        s = nf90_put_var(nc, idLon,  totglamt); call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idLat,  totgphit); call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idLev,     gdept); call handle_err1(s,counter,fileNetCDF)

        s = nf90_put_var(nc, idTim,    julian); call handle_err1(s,counter,fileNetCDF)

        s = nf90_put_var(nc, idB,      tottrb); call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idN,      tottrn); call handle_err1(s,counter,fileNetCDF)

        s =nf90_close(nc)


       END SUBROUTINE write_restart_ben

       !****************************************************************************
       !****************************************************************************
       !****************************************************************************

       SUBROUTINE WRITE_AVE_ben(fileNetCDF, datefrom, dateTo,M)
       USE netcdf
       USE BIOBEN_mem
       USE myalloc
       IMPLICIT NONE

       CHARACTER*(*) fileNetCDF
       character(LEN=17) datefrom, dateTo
       real(4) M(jpiglo, jpjglo, jpk_ben)

       character(LEN=20) VAR
       integer istart,iend

       integer s, nc, counter
       integer timid, depid, yid, xid
       integer idvartime,idgdept,idphit,idlamt,idVAR
       real lat_actual_range(2), lon_actual_range(2), depth_actual_range(2)
         lon_actual_range=(/-9.25  , 36.0   /)
         lat_actual_range=(/30.5   , 44.5   /)
       depth_actual_range=(/ 4.9991,4450.068/)

       istart=index(fileNetCDF,'00.')+3
       iend  =index(fileNetCDF,'.nc')-1
       VAR        =fileNetCDF(istart:iend)


        s = nf90_create(fileNetCDF, NF90_CLOBBER, nc)
        ! *********** GLOBAL ********************
        s = nf90_put_att(nc, nf90_global, 'Convenctions'   ,'COARDS')
        s = nf90_put_att(nc, nf90_global, 'DateStart'     , datefrom)
        s = nf90_put_att(nc, nf90_global, 'Date__End'     ,   dateTo)

        ! *********** DIMENSIONS ****************
        s= nf90_def_dim(nc,'lon'           , jpiglo,  xid)
        s= nf90_def_dim(nc,'lat'           , jpjglo,  yid)
        s= nf90_def_dim(nc,'depth'         , jpk_ben   ,depid)
        s= nf90_def_dim(nc,'time'  , NF90_UNLIMITED,timid)

        ! ********** VARIABLES *****************
        s = nf90_def_var(nc,'time',         nf90_double,(/timid/),       idvartime)
        s = nf90_def_var(nc,'depth',        nf90_float, (/depid/),         idgdept)
        s = nf90_def_var(nc,'lat'   ,       nf90_float, (/yid/),            idphit)
        s = nf90_def_var(nc,'lon'   ,       nf90_float, (/xid/),            idlamt)

       s = nf90_def_var(nc,trim(VAR) ,        nf90_float, (/xid,yid,depid,timid/),  idVAR)

        s = nf90_put_att(nc,idgdept,'units'        ,'meter')
        s = nf90_put_att(nc,idgdept,'positive'     ,'down')
        s = nf90_put_att(nc,idgdept,'actual_range' ,depth_actual_range)

        s = nf90_put_att(nc,idphit, 'units'        ,'degrees_north')
        s = nf90_put_att(nc,idphit, 'long_name'    ,'Latitude')
        s = nf90_put_att(nc,idphit,'actual_range' ,lat_actual_range)

        s = nf90_put_att(nc,idlamt, 'units'        ,'degrees_east')
        s = nf90_put_att(nc,idlamt, 'long_name'    ,'Longitude')
        s = nf90_put_att(nc,idlamt,'actual_range' ,lon_actual_range)

        s = nf90_put_att(nc,idVAR, 'long_name'    ,VAR)
        s = nf90_put_att(nc,idVAR, 'missing_value'    ,1.e+20)

        s =nf90_enddef(nc)

        counter=0

        s = nf90_put_var(nc, idlamt,   REAL(totglamt(:,jpjglo),4) ); call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idphit,   REAL(totgphit(jpiglo,:),4) ); call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idgdept,  REAL(   gdept,     4) ); call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idVAR  ,  M )                    ; call handle_err1(s,counter,fileNetCDF)


        s =nf90_close(nc)

       END SUBROUTINE WRITE_AVE_ben

       SUBROUTINE WRITE_AVE_2D_ben(fileNetCDF, datefrom, dateTo,M)
       USE netcdf
       USE BIOBEN_mem
       USE myalloc
       IMPLICIT NONE

       CHARACTER*(*) fileNetCDF
       character(LEN=17) datefrom, dateTo
       real(4) M(jpiglo, jpjglo)

       character(LEN=20) VAR
       integer istart,iend

       integer s, nc, counter
       integer timid, yid, xid
       integer idvartime,idphit,idlamt,idVAR
       real lat_actual_range(2), lon_actual_range(2)
         lon_actual_range=(/-9.25  , 36.0   /)
         lat_actual_range=(/30.5   , 44.5   /)


       istart=index(fileNetCDF,'00.')+3
       iend  =index(fileNetCDF,'.nc')-1
       VAR        =fileNetCDF(istart:iend)


        s = nf90_create(fileNetCDF, NF90_CLOBBER, nc)
        ! *********** GLOBAL ********************
        s = nf90_put_att(nc, nf90_global, 'Convenctions'   ,'COARDS')
        s = nf90_put_att(nc, nf90_global, 'DateStart'     , datefrom)
        s = nf90_put_att(nc, nf90_global, 'Date__End'     ,   dateTo)

        ! *********** DIMENSIONS ****************
        s= nf90_def_dim(nc,'lon'           , jpiglo,  xid)
        s= nf90_def_dim(nc,'lat'           , jpjglo,  yid)
        s= nf90_def_dim(nc,'time'  , NF90_UNLIMITED,timid)

        ! ********** VARIABLES *****************
        s = nf90_def_var(nc,'time',         nf90_double,(/timid/),       idvartime)
        s = nf90_def_var(nc,'lat'   ,       nf90_float, (/yid/),            idphit)
        s = nf90_def_var(nc,'lon'   ,       nf90_float, (/xid/),            idlamt)

       s = nf90_def_var(nc,trim(VAR) ,        nf90_float, (/xid,yid,timid/),  idVAR)

        s = nf90_put_att(nc,idphit, 'units'        ,'degrees_north')
        s = nf90_put_att(nc,idphit, 'long_name'    ,'Latitude')
        s = nf90_put_att(nc,idphit,'actual_range' ,lat_actual_range)

        s = nf90_put_att(nc,idlamt, 'units'        ,'degrees_east')
        s = nf90_put_att(nc,idlamt, 'long_name'    ,'Longitude')
        s = nf90_put_att(nc,idlamt,'actual_range' ,lon_actual_range)

        s = nf90_put_att(nc,idVAR, 'long_name'    ,VAR)
        s = nf90_put_att(nc,idVAR, 'missing_value' ,1.e+20)

        s =nf90_enddef(nc)

        counter=0

        s = nf90_put_var(nc, idlamt,   REAL(totglamt(:,jpjglo),4) ); call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idphit,   REAL(totgphit(jpiglo,:),4) ); call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idVAR  ,  M )                    ; call handle_err1(s,counter,fileNetCDF)


        s =nf90_close(nc)

       END SUBROUTINE WRITE_AVE_2D_ben

!****************************************************************************
       !****************************************************************************
       !****************************************************************************

       SUBROUTINE WRITE_AVE_BKP_ben(fileNetCDF, datefrom, dateTo,M, ave_counter)
       USE netcdf
       USE BIOBEN_mem
       USE myalloc
       IMPLICIT NONE

       CHARACTER*(*) fileNetCDF
       character(LEN=17) datefrom, dateTo
       real(8) M(jpiglo, jpjglo, jpk_ben)
       integer ave_counter

       !local
       character(LEN=20) VAR
       integer istart,iend

       integer s, nc, counter
       integer timid, depid, yid, xid
       integer idvartime,idgdept,idphit,idlamt,idVAR


       istart=index(fileNetCDF,'00.')+3
       iend  =index(fileNetCDF,'.nc')-1
       VAR        =fileNetCDF(istart:iend)


        s = nf90_create(fileNetCDF, NF90_CLOBBER, nc)
        ! *********** GLOBAL ********************
        s = nf90_put_att(nc, nf90_global, 'Convenctions'  ,    'COARDS')
        s = nf90_put_att(nc, nf90_global, 'DateStart'     ,    datefrom)
        s = nf90_put_att(nc, nf90_global, 'Date__End'     ,      dateTo)
        s = nf90_put_att(nc, nf90_global, 'ave_counter'   , ave_counter)

        ! *********** DIMENSIONS ****************
        s= nf90_def_dim(nc,'lon'           , jpiglo,  xid)
        s= nf90_def_dim(nc,'lat'           , jpjglo,  yid)
        s= nf90_def_dim(nc,'depth'         , jpk_ben   ,depid)
        s= nf90_def_dim(nc,'time'  , NF90_UNLIMITED,timid)

        ! ********** VARIABLES *****************
        s = nf90_def_var(nc,'time',         nf90_double,(/timid/),       idvartime)
        s = nf90_def_var(nc,'depth',        nf90_float, (/depid/),         idgdept)
        s = nf90_def_var(nc,'lat'   ,       nf90_float, (/yid/),            idphit)
        s = nf90_def_var(nc,'lon'   ,       nf90_float, (/xid/),            idlamt)

       s = nf90_def_var(nc,trim(VAR) ,           nf90_double,(/xid,yid,depid,timid/),  idVAR)

        s = nf90_put_att(nc,idgdept,'units'        ,'m')
        s = nf90_put_att(nc,idgdept,'positive'     ,'down')

        s = nf90_put_att(nc,idphit, 'units'        ,'degrees_north')
        s = nf90_put_att(nc,idphit, 'long_name'    ,'Latitude')

        s = nf90_put_att(nc,idlamt, 'units'        ,'degrees_east')
        s = nf90_put_att(nc,idlamt, 'long_name'    ,'Longitude')


        s = nf90_put_att(nc,idVAR, 'long_name'        ,VAR)
        s = nf90_put_att(nc,idVAR, 'missing_value'    ,1.e+20)

        s =nf90_enddef(nc)

        counter=0


        s = nf90_put_var(nc, idlamt,   REAL(totglamt(:,jpjglo),4) ); call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idphit,   REAL(totgphit(jpiglo,:),4) ); call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idgdept,  REAL(   gdept,          4) ); call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idVAR  ,                          M  ); call handle_err1(s,counter,fileNetCDF)


        s =nf90_close(nc)

       END SUBROUTINE WRITE_AVE_BKP_ben

      SUBROUTINE WRITE_AVE_2d_BKP_ben(fileNetCDF, datefrom, dateTo,M, ave_counter)
       USE netcdf
       USE BIOBEN_mem
       USE myalloc
       IMPLICIT NONE

       CHARACTER*(*) fileNetCDF
       character(LEN=17) datefrom, dateTo
       real(8) M(jpiglo, jpjglo)
       integer ave_counter

       !local
       character(LEN=20) VAR
       integer istart,iend

       integer s, nc, counter
       integer timid, yid, xid
       integer idvartime,idgdept,idphit,idlamt,idVAR


       istart=index(fileNetCDF,'00.')+3
       iend  =index(fileNetCDF,'.nc')-1
       VAR        =fileNetCDF(istart:iend)


        s = nf90_create(fileNetCDF, NF90_CLOBBER, nc)
        ! *********** GLOBAL ********************
        s = nf90_put_att(nc, nf90_global, 'Convenctions'  ,    'COARDS')
        s = nf90_put_att(nc, nf90_global, 'DateStart'     ,    datefrom)
        s = nf90_put_att(nc, nf90_global, 'Date__End'     ,      dateTo)
        s = nf90_put_att(nc, nf90_global, 'ave_counter'   , ave_counter)

        ! *********** DIMENSIONS ****************
        s= nf90_def_dim(nc,'lon'           , jpiglo,  xid)
        s= nf90_def_dim(nc,'lat'           , jpjglo,  yid)
        s= nf90_def_dim(nc,'time'  , NF90_UNLIMITED,timid)

        ! ********** VARIABLES *****************
        s = nf90_def_var(nc,'time',         nf90_double,(/timid/),       idvartime)
        s = nf90_def_var(nc,'lat'   ,       nf90_float, (/yid/),            idphit)
        s = nf90_def_var(nc,'lon'   ,       nf90_float, (/xid/),            idlamt)

       s = nf90_def_var(nc,trim(VAR) ,           nf90_double,(/xid,yid,timid/),  idVAR)

        s = nf90_put_att(nc,idphit, 'units'        ,'degrees_north')
        s = nf90_put_att(nc,idphit, 'long_name'    ,'Latitude')

        s = nf90_put_att(nc,idlamt, 'units'        ,'degrees_east')
        s = nf90_put_att(nc,idlamt, 'long_name'    ,'Longitude')


        s = nf90_put_att(nc,idVAR, 'long_name'        ,VAR)
        s = nf90_put_att(nc,idVAR, 'missing_value'    ,1.e+20)

        s =nf90_enddef(nc)

        counter=0


        s = nf90_put_var(nc, idlamt,   REAL(totglamt(:,jpjglo),4) ); call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idphit,   REAL(totgphit(jpiglo,:),4) ); call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idVAR  ,                          M  ); call handle_err1(s,counter,fileNetCDF)


        s =nf90_close(nc)

       END SUBROUTINE WRITE_AVE_2d_BKP_ben



       !****************************************************************************

       !****************************************************************************
       !****************************************************************************


        SUBROUTINE get_att_char_ben(fileNetCDF,attname,attvalue)
        use netcdf
        IMPLICIT NONE
        character fileNetCDF*(*)
        character attname*(*)
        character attvalue*(*)

        integer stat, ncid
        integer counter
        counter=0
        stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)  ;       call handle_err1(stat, counter,FileNetCDF)
        stat= nf90_get_att(ncid, nf90_global, attname, attvalue); call handle_err1(stat, counter,FileNetCDF)
        stat = nf90_close(ncid);

        END SUBROUTINE get_att_char_ben



        !****************************************************************************
        !****************************************************************************

        SUBROUTINE get_att_int_ben(fileNetCDF,attname,attvalue)
        use netcdf
        IMPLICIT NONE
        character fileNetCDF*(*)
        character attname*(*)
        integer attvalue

        integer stat, ncid
        integer counter
        counter=0
        stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)  ;       call handle_err1(stat, counter,FileNetCDF)
        stat= nf90_get_att(ncid, nf90_global, attname, attvalue); call handle_err1(stat, counter,FileNetCDF)
        stat = nf90_close(ncid);

        END SUBROUTINE get_att_int_ben


      !****************************************************************************
        subroutine handle_err1_ben(status,mycount, fileNetCDF)
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

        end subroutine handle_err1_ben


      ! **************************************************************************
        subroutine handle_err2_ben(status,fileNetCDF,varname)
        USE netcdf
        integer status
        character fileNetCDF*(*) ,varname*(*)

        if(status .ne. nf90_NoErr)  then
           write(*,*) 'ERROR in Var = ', varname, ' file :', fileNetCDF
        endif

        end subroutine handle_err2_ben

       ! ***************************************************************************



