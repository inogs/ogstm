        SUBROUTINE CREATEMISFIT(SATFILE,VARFILE,OPTION, ISLOG, MISFIT_FILE)

        !createMisfit $DA/BILINEAR/sat_file_8.nc $PWD/MODEL_2D_chl.nc $VAR2D_FILE $OPT $ISLOG
        ! OPT can assume values 1,2,3
        ! ISLOG can assume values 0 or 1

          use DA_MEM
          use netcdf

          implicit none


          integer i,j
          character*1024  SATFILE, VARFILE, MISFIT_FILE
          character       OPTION, ISLOGCHAR
          character description*1024
          logical B, ISLOG
          real fillValue, ERRSATuniformVALUE



           write(*,*) 'SATFILE =', trim(SATFILE)
           write(*,*) 'VARFILE =', trim(VARFILE)
           if (ISLOG) then
                 write(*,*) 'Logarithms of concentrations will be taken in account'
                 ERRSATuniformVALUE = 0.3
            else
                 write(*,*) 'Concentrations will be taken in account'
                 ERRSATuniformVALUE = 0.73
           END if


          SELECT CASE(OPTION)
            CASE('1')
                write(description,'(A,ES11.4)') 'ERRsat = uniform value over the Mediterranean Sea', ERRSATuniformVALUE
            CASE('2')
                description ='ERRsat read from file: '//trim(VARFILE)
            CASE('3')
                description ='ERRsat = LOG(.01)+ .3*LOG(CHLsat)'
          END SELECT

           write(*,*) trim(description)





          call readNetCDF_2dvar(SATFILE,trim(satvarname)   ,jpiglo,jpjglo,  CHLsat  ) ! testcase value around 0.25
          call readNetCDF_2dvar(VARFILE,'variance',jpiglo,jpjglo,  VAR2D   ) !                around 0.0005


          fillValue = 1.0e+20

          do i=1,jpiglo
          do j=1,jpjglo


            if ( isnan2(CHLsat(i,j)).or.(CHLsat(i,j).eq.fillValue)) then
               ERRsat(i,j) = fillValue
               MISFIT(i,j) = fillValue
              cMISFIT(i,j) = fillValue
            else
              !****************************
        !       if ( ERRsat(i,j).gt.sqrt(VAR2D(i,j)) ) ERRsat(i,j) =  sqrt(VAR2D(i,j))

              SELECT CASE (OPTION)
                CASE('1')
                    ERRsat(i,j) = ERRSATuniformVALUE
                CASE('2')
                    ERRsat(i,j) = sqrt(VAR2D(i,j))
                    if (ERRsat(i,j).gt.ERRSATuniformVALUE) ERRsat(i,j) = ERRSATuniformVALUE
                CASE('3')
                    ERRsat(i,j)=  LOG(.01)+ .3*LOG(CHLsat(i,j))
              END SELECT

              !****************************

               if (CHL_SUP(i,j).lt.fillValue) then
                !***************************************
                  MISFIT(i,j) = LOG(CHLsat(i,j)) - LOG(CHL_SUP(i,j));
                 cMISFIT(i,j) =     CHLsat(i,j)      - CHL_SUP(i,j);
                !***************************************
               else
                  MISFIT(i,j) = fillValue
                 cMISFIT(i,j) = fillValue
               endif

            endif

          end do
          end do


        ! NaN number of VAR2D must result in a NaN in MISFIT
          if (OPTION.eq.'2') then
             do i=1,jpiglo
             do j=1,jpjglo
               if (VAR2D(i,j).eq.fillValue) then
                    MISFIT(i,j) = fillValue
                   cMISFIT(i,j) = fillValue
               endif
             enddo
             enddo
          endif


          if (.not.ISLOG) MISFIT=cMISFIT



        if (ISLOG) then
           B= writeMISFIT_LOG(MISFIT_FILE)
        else
           B= writeMISFIT(MISFIT_FILE)
        endif

        write(*,*) trim(MISFIT_FILE)//' created'

        CONTAINS

        ! ***********************************************************************
        LOGICAL FUNCTION isnan2(A)
        implicit none
        real, INTENT(IN) :: A
        if ( A.eq.A ) then
          isnan2 = .FALSE.
          else
          isnan2 = .TRUE.
        end if
        END FUNCTION isnan2


        ! ***********************************************************************
        ! ***********************************************************************

        LOGICAL FUNCTION writeMISFIT_LOG(fileNetCDF)
        use myalloc
        use netcdf
        implicit none

        character, INTENT(IN) :: fileNetCDF*(*)


        ! local
        integer ncid, s, counter
        integer IDtime, IDdepth, IDlon, IDlat
        integer ID_time,ID_depth,ID_lon,ID_lat, ID_misf,ID_err, ID_Nolog
        real fillValue, depth
        real time

        fillValue=1.0e+20
        time     =1.
        depth    =1.47210180759
        counter=0

        s = nf90_create(path = fileNetCDF, cmode= NF90_CLOBBER , ncid = ncid)
        call handle_err1(s, counter,FileNetCDF)



        s = nf90_def_dim(ncid, 'time' , 1, IDtime)  ; call handle_err1(s, counter,FileNetCDF)
        s = nf90_def_dim(ncid, 'depht', 1, IDdepth) ; call handle_err1(s, counter,FileNetCDF)
        s = nf90_def_dim(ncid, 'lon'  , jpiglo, IDlon)  ; call handle_err1(s, counter,FileNetCDF)
        s = nf90_def_dim(ncid, 'lat'  , jpjglo, IDlat)  ; call handle_err1(s, counter,FileNetCDF)

        s = nf90_def_var(ncid, 'misfchl'    , nf90_float, (/IDlon,IDlat/), ID_misf ); call handle_err1(s, counter,FileNetCDF)
        s = nf90_def_var(ncid, 'errchl'     , nf90_float, (/IDlon,IDlat/), ID_err  ); call handle_err1(s, counter,FileNetCDF)
        s = nf90_def_var(ncid, 'misfit_conc', nf90_float, (/IDlon,IDlat/), ID_Nolog); call handle_err1(s, counter,FileNetCDF)
        s = nf90_def_var(ncid, 'depth'      , nf90_float, (/IDdepth/)    , ID_depth); call handle_err1(s, counter,FileNetCDF)
        s = nf90_def_var(ncid, 'time'       , nf90_float, (/IDtime/)     , ID_time ); call handle_err1(s, counter,FileNetCDF)
        s = nf90_def_var(ncid, 'lon'        , nf90_float, (/IDlon/)      , ID_lon  ); call handle_err1(s, counter,FileNetCDF)
        s = nf90_def_var(ncid, 'lat'        , nf90_float, (/IDlat/)      , ID_lat  ); call handle_err1(s, counter,FileNetCDF)



        s = nf90_put_att(ncid,ID_lon,   '_FillValue',fillValue); call handle_err1(s, counter,FileNetCDF)
        s = nf90_put_att(ncid,ID_lat,   '_FillValue',fillValue); call handle_err1(s, counter,FileNetCDF)

        s = nf90_put_att(ncid,ID_err,   '_FillValue',fillValue); call handle_err1(s, counter,FileNetCDF)
        s = nf90_put_att(ncid,ID_err,   'info'    ,description); call handle_err1(s, counter,FileNetCDF)

        s = nf90_put_att(ncid,ID_misf,  '_FillValue',fillValue); call handle_err1(s, counter,FileNetCDF)
        s = nf90_put_att(ncid,ID_misf,  'type','logarithm of chl concentrations'); call handle_err1(s, counter,FileNetCDF)

        s = nf90_put_att(ncid,ID_Nolog, '_FillValue',fillValue); call handle_err1(s, counter,FileNetCDF)
        s = nf90_put_att(ncid,ID_Nolog, 'type','chl concentrations'); call handle_err1(s, counter,FileNetCDF)
        s = nf90_put_att(ncid,ID_depth, '_FillValue',fillValue); call handle_err1(s, counter,FileNetCDF)
        s = nf90_put_att(ncid,ID_time,  '_FillValue',fillValue); call handle_err1(s, counter,FileNetCDF)


        s = nf90_enddef(ncid)
        call handle_err1(s, counter,FileNetCDF)


        s = nf90_put_var(ncid, ID_lon   , REAL(totglamt(:,jpjglo),4)    ); call handle_err1(s, counter,FileNetCDF)
        s = nf90_put_var(ncid, ID_lat   , REAL(totgphit(jpiglo,:),4)    ); call handle_err1(s, counter,FileNetCDF)
        s = nf90_put_var(ncid, ID_misf  , MISFIT ); call handle_err1(s, counter,FileNetCDF)
        s = nf90_put_var(ncid, ID_err   , ERRsat ); call handle_err1(s, counter,FileNetCDF)
        s = nf90_put_var(ncid, ID_NoLog , cMISFIT); call handle_err1(s, counter,FileNetCDF)
        s = nf90_put_var(ncid, ID_time  , time   ); call handle_err1(s, counter,FileNetCDF)
        s = nf90_put_var(ncid, ID_depth , depth  ); call handle_err1(s, counter,FileNetCDF)


        s=nf90_close(ncid)
        writeMISFIT_LOG=.true.
        END FUNCTION writeMISFIT_LOG


        ! ***********************************************************************
        ! ***********************************************************************

        LOGICAL FUNCTION writeMISFIT(fileNetCDF)
        use netcdf
        use myalloc
        implicit none

        character, INTENT(IN) :: fileNetCDF*(*)


        ! local
        integer ncid, s, counter
        integer IDtime, IDdepth, IDlon, IDlat
        integer ID_time,ID_depth,ID_lon,ID_lat, ID_misf,ID_err, ID_Nolog
        real fillValue, depth
        real time

        fillValue=1.0e+20
        time     =1.
        depth    =1.47210180759
        counter = 0

        s = nf90_create(path = fileNetCDF, cmode= NF90_CLOBBER , ncid = ncid)
        call handle_err1(s, counter,FileNetCDF)


        s = nf90_def_dim(ncid, 'time' , 1, IDtime)  ; call handle_err1(s, counter,FileNetCDF)
        s = nf90_def_dim(ncid, 'depht', 1, IDdepth) ; call handle_err1(s, counter,FileNetCDF)
        s = nf90_def_dim(ncid, 'lon'  , jpiglo, IDlon)  ; call handle_err1(s, counter,FileNetCDF)
        s = nf90_def_dim(ncid, 'lat'  , jpjglo, IDlat)  ; call handle_err1(s, counter,FileNetCDF)

        s = nf90_def_var(ncid, 'misfchl'    , nf90_float, (/IDlon,IDlat/), ID_misf ); call handle_err1(s, counter,FileNetCDF)
        s = nf90_def_var(ncid, 'errchl'     , nf90_float, (/IDlon,IDlat/), ID_err  ); call handle_err1(s, counter,FileNetCDF)
        s = nf90_def_var(ncid, 'depth'      , nf90_float, (/IDdepth/)    , ID_depth); call handle_err1(s, counter,FileNetCDF)
        s = nf90_def_var(ncid, 'time'       , nf90_float, (/IDtime/)     , ID_time ); call handle_err1(s, counter,FileNetCDF)
        s = nf90_def_var(ncid, 'lon'        , nf90_float, (/IDlon/)      , ID_lon  ); call handle_err1(s, counter,FileNetCDF)
        s = nf90_def_var(ncid, 'lat'        , nf90_float, (/IDlat/)      , ID_lat  ); call handle_err1(s, counter,FileNetCDF)



        s = nf90_put_att(ncid,ID_lon,   '_FillValue',fillValue); call handle_err1(s, counter,FileNetCDF)
        s = nf90_put_att(ncid,ID_lat,   '_FillValue',fillValue); call handle_err1(s, counter,FileNetCDF)

        s = nf90_put_att(ncid,ID_err,   '_FillValue',fillValue); call handle_err1(s, counter,FileNetCDF)
        s = nf90_put_att(ncid,ID_err,   'info'    ,description); call handle_err1(s, counter,FileNetCDF)

        s = nf90_put_att(ncid,ID_misf,  '_FillValue',fillValue); call handle_err1(s, counter,FileNetCDF)
        s = nf90_put_att(ncid,ID_misf,  'type','chl concentrations'); call handle_err1(s, counter,FileNetCDF)

        s = nf90_put_att(ncid,ID_depth, '_FillValue',fillValue); call handle_err1(s, counter,FileNetCDF)
        s = nf90_put_att(ncid,ID_time,  '_FillValue',fillValue); call handle_err1(s, counter,FileNetCDF)


        s = nf90_enddef(ncid)
        call handle_err1(s, counter,FileNetCDF)

        s = nf90_put_var(ncid, ID_lon   , REAL(totglamt(jpjglo,:),4)    ); call handle_err1(s, counter,FileNetCDF)
        s = nf90_put_var(ncid, ID_lat   , REAL(totgphit(:,jpiglo),4)    ); call handle_err1(s, counter,FileNetCDF)
        s = nf90_put_var(ncid, ID_misf  , MISFIT ); call handle_err1(s, counter,FileNetCDF)
        s = nf90_put_var(ncid, ID_err   , ERRsat ); call handle_err1(s, counter,FileNetCDF)
        s = nf90_put_var(ncid, ID_time  , time   ); call handle_err1(s, counter,FileNetCDF)
        s = nf90_put_var(ncid, ID_depth , depth  ); call handle_err1(s, counter,FileNetCDF)


        s=nf90_close(ncid)
        writeMISFIT=.true.
        END FUNCTION writeMISFIT

        END SUBROUTINE CREATEMISFIT
