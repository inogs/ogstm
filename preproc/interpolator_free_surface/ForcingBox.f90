    LOGICAL FUNCTION PREPARE_BOX_FILE_U(fileNetCDF)
        use netcdf
        implicit none

        character fileNetCDF*(*)
        integer s, nc,G, counter
        integer timid, depid, yid, xid
        integer idvartime, idgdept, idphit, idlamt, idU

        G     = nf90_global

        s = nf90_create(fileNetCDF, NF90_CLOBBER, nc)
        ! *********** GLOBAL ********************
        s = nf90_put_att(nc, G, 'Convenctions'     , 'OPA')


        ! *********** DIMENSIONS ****************
        s= nf90_def_dim(nc,'x'           , im,  xid)
        s= nf90_def_dim(nc,'y'           , jm,  yid)
        s= nf90_def_dim(nc,'deptht'      , km,depid)
        s= nf90_def_dim(nc,'time_counter', NF90_UNLIMITED,timid)


        ! ********** VARIABLES *****************
        s = nf90_def_var(nc,'time_counter',  nf90_double,(/timid/),           idvartime)
        s = nf90_def_var(nc,'deptht',        nf90_float, (/depid/),             idgdept)
        s = nf90_def_var(nc,'nav_lat',       nf90_float, (/xid,yid/),            idphit)
        s = nf90_def_var(nc,'nav_lon',       nf90_float, (/xid,yid/),            idlamt)
        s = nf90_def_var(nc,'vozocrtx' ,     nf90_float, (/xid,yid,depid,timid/),   idU)



        s = nf90_put_att(nc,idvartime,'units', UnitsTime )
        s = nf90_put_att(nc,idvartime,'calendar', 'gregorian' )
        s = nf90_put_att(nc,idgdept,'units'        ,'m')
        s = nf90_put_att(nc,idgdept,'positive'     ,'down')
        s = nf90_put_att(nc,idphit, 'units'        ,'degrees_north')
        s = nf90_put_att(nc,idphit, 'long_name'    ,'Latitude')
        s = nf90_put_att(nc,idlamt, 'units'        ,'degrees_east')
        s = nf90_put_att(nc,idlamt, 'long_name'    ,'Longitude')

        s = nf90_put_att(nc,idU   , 'long_name'    ,'Zonal Current')
        s = nf90_put_att(nc,idU   , 'units'         ,'m/s')
        s = nf90_put_att(nc,idU   , 'missing_value',1.e+20)

        s =nf90_enddef(nc)

        counter=0

        s = nf90_put_var(nc, idlamt,       X  );call handle_err1(s,counter)
        s = nf90_put_var(nc, idphit,       Y  );call handle_err1(s,counter)
        s = nf90_put_var(nc, idgdept,      Z  );call handle_err1(s,counter)
        s = nf90_put_var(nc, idvartime, TheTime  );
        call handle_err1(s,counter)

        s= nf90_close(nc)
        PREPARE_BOX_FILE_U = .true.
    END FUNCTION PREPARE_BOX_FILE_U



    ! **************************************************************************
    ! **************************************************************************

    LOGICAL FUNCTION PREPARE_BOX_FILE_V(fileNetCDF)
        use netcdf
        implicit none

        character fileNetCDF*(*)
        integer s, nc,G, counter
        integer timid, depid, yid, xid
        integer idvartime, idgdept, idphit, idlamt, idV

        G     = nf90_global

        s = nf90_create(fileNetCDF, NF90_CLOBBER, nc)
        ! *********** GLOBAL ********************
        s = nf90_put_att(nc, G, 'Convenctions'     , 'OPA')

        ! *********** DIMENSIONS ****************
        s= nf90_def_dim(nc,'x'           , im,  xid)
        s= nf90_def_dim(nc,'y'           , jm,  yid)
        s= nf90_def_dim(nc,'deptht'      , km,depid)
        s= nf90_def_dim(nc,'time_counter', NF90_UNLIMITED,timid)


        ! ********** VARIABLES *****************
        s = nf90_def_var(nc,'time_counter',  nf90_double,(/timid/),           idvartime)
        s = nf90_def_var(nc,'deptht',        nf90_float, (/depid/),             idgdept)
        s = nf90_def_var(nc,'nav_lat',       nf90_float, (/xid,yid/),            idphit)
        s = nf90_def_var(nc,'nav_lon',       nf90_float, (/xid,yid/),            idlamt)
        s = nf90_def_var(nc,'vomecrty' ,     nf90_float, (/xid,yid,depid,timid/),   idV)



        s = nf90_put_att(nc,idvartime,'units', UnitsTime )
        s = nf90_put_att(nc,idgdept,'units'        ,'m')
        s = nf90_put_att(nc,idgdept,'positive'     ,'down')
        s = nf90_put_att(nc,idphit, 'units'        ,'degrees_north')
        s = nf90_put_att(nc,idphit, 'long_name'    ,'Latitude')
        s = nf90_put_att(nc,idlamt, 'units'        ,'degrees_east')
        s = nf90_put_att(nc,idlamt, 'long_name'    ,'Longitude')

        s = nf90_put_att(nc,idV   , 'long_name'    ,'Meridional Current')
        s = nf90_put_att(nc,idV   , 'units'        ,'m/s')
        s = nf90_put_att(nc,idV   , 'missing_value',1.e+20)

        s =nf90_enddef(nc)

        counter=0

        s = nf90_put_var(nc, idlamt,       X  );call handle_err1(s,counter)
        s = nf90_put_var(nc, idphit,       Y  );call handle_err1(s,counter)
        s = nf90_put_var(nc, idgdept,      Z  );call handle_err1(s,counter)
        s = nf90_put_var(nc, idvartime, TheTime  );
        call handle_err1(s,counter)

        s= nf90_close(nc)
        PREPARE_BOX_FILE_V = .true.
    END FUNCTION PREPARE_BOX_FILE_V

    ! **************************************************************************
    ! **************************************************************************

    LOGICAL FUNCTION PREPARE_BOX_FILE_W(fileNetCDF)
        use netcdf
        implicit none

        character fileNetCDF*(*)
        integer s, nc,G, counter
        integer timid, depid, yid, xid
        integer idvartime, idgdept, idphit, idlamt, idW, idEddy

        G     = nf90_global

        s = nf90_create(fileNetCDF, NF90_CLOBBER, nc)
        ! *********** GLOBAL ********************
        s = nf90_put_att(nc, G, 'Convenctions'     , 'OPA')

        ! *********** DIMENSIONS ****************
        s= nf90_def_dim(nc,'x'           , im,  xid)
        s= nf90_def_dim(nc,'y'           , jm,  yid)
        s= nf90_def_dim(nc,'depthw'      , km,depid)
        s= nf90_def_dim(nc,'time_counter', NF90_UNLIMITED,timid)


        ! ********** VARIABLES *****************
        s = nf90_def_var(nc,'time_counter',  nf90_double,(/timid/),           idvartime)
        s = nf90_def_var(nc,'depthw',        nf90_float, (/depid/),             idgdept)
        s = nf90_def_var(nc,'nav_lat',       nf90_float, (/xid,yid/),            idphit)
        s = nf90_def_var(nc,'nav_lon',       nf90_float, (/xid,yid/),            idlamt)
        s = nf90_def_var(nc,'vovecrtz' ,     nf90_float, (/xid,yid,depid,timid/),   idW)
        s = nf90_def_var(nc,'votkeavt' ,     nf90_float, (/xid,yid,depid,timid/),idEddy)


        s = nf90_put_att(nc,idvartime,'units', UnitsTime )
        s = nf90_put_att(nc,idgdept,'units'        ,'m')
        s = nf90_put_att(nc,idgdept,'positive'     ,'down')
        s = nf90_put_att(nc,idphit, 'units'        ,'degrees_north')
        s = nf90_put_att(nc,idphit, 'long_name'    ,'Latitude')
        s = nf90_put_att(nc,idlamt, 'units'        ,'degrees_east')
        s = nf90_put_att(nc,idlamt, 'long_name'    ,'Longitude')

        s = nf90_put_att(nc,idW   , 'long_name'    ,'Vertical Velocity')
        s = nf90_put_att(nc,idW   , 'units'        ,'m/s')
        s = nf90_put_att(nc,idW   , 'missing_value',1.e+20)

        s = nf90_put_att(nc,idEddy, 'long_name'    ,'Vertical Eddy Diffusivity')
        s = nf90_put_att(nc,idEddy, 'units'        ,'m2/s')
        s = nf90_put_att(nc,idEddy, 'missing_value',1.e+20)

        s =nf90_enddef(nc)

        counter=0

        s = nf90_put_var(nc, idlamt,       X  );call handle_err1(s,counter)
        s = nf90_put_var(nc, idphit,       Y  );call handle_err1(s,counter)
        s = nf90_put_var(nc, idgdept,      Z  );call handle_err1(s,counter)
        s = nf90_put_var(nc, idvartime, TheTime  );
        call handle_err1(s,counter)

        s= nf90_close(nc)
        PREPARE_BOX_FILE_W = .true.
    END FUNCTION PREPARE_BOX_FILE_W

    ! **************************************************************************
    ! **************************************************************************

    LOGICAL FUNCTION PREPARE_BOX_FILE_T(fileNetCDF)
        use netcdf
        implicit none

        character fileNetCDF*(*)
        integer s, nc,G, counter
        integer timid, depid, yid, xid
        integer idvartime, idgdept, idphit, idlamt, idT, idS, idR, idW,idH,idE

        G     = nf90_global

        s = nf90_create(fileNetCDF, NF90_CLOBBER, nc)
        ! *********** GLOBAL ********************
        s = nf90_put_att(nc, G, 'Convenctions'     , 'OPA')

        ! *********** DIMENSIONS ****************
        s= nf90_def_dim(nc,'x'           , im,  xid)
        s= nf90_def_dim(nc,'y'           , jm,  yid)
        s= nf90_def_dim(nc,'deptht'      , km,depid)
        s= nf90_def_dim(nc,'time_counter', NF90_UNLIMITED,timid)


        ! ********** VARIABLES *****************
        s = nf90_def_var(nc,'time_counter',  nf90_double,(/timid/),           idvartime)
        s = nf90_def_var(nc,'deptht',        nf90_float, (/depid/),             idgdept)
        s = nf90_def_var(nc,'nav_lat',       nf90_float, (/xid,yid/),            idphit)
        s = nf90_def_var(nc,'nav_lon',       nf90_float, (/xid,yid/),            idlamt)
        s = nf90_def_var(nc,'vosaline' ,     nf90_float, (/xid,yid,depid,timid/),   idS)
        s = nf90_def_var(nc,'votemper' ,     nf90_float, (/xid,yid,depid,timid/),   idT)
        s = nf90_def_var(nc,'soshfldo' ,     nf90_float, (/xid,yid,      timid/),   idR)
        s = nf90_def_var(nc,'sowindsp' ,     nf90_float, (/xid,yid,      timid/),   idW)
        s = nf90_def_var(nc,'sossheig' ,     nf90_float, (/xid,yid,      timid/),   idH)
        s = nf90_def_var(nc,'sowaflcd' ,     nf90_float, (/xid,yid,      timid/),   idE)

        s = nf90_put_att(nc,idvartime,'units', UnitsTime )

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

        s = nf90_put_att(nc,idR   , 'long_name'    ,'Short_Wave_Radiation')
        s = nf90_put_att(nc,idR   , 'units'        ,'W/m2')
        s = nf90_put_att(nc,idR   , 'missing_value',1.e+20)

        s = nf90_put_att(nc,idW   , 'long_name'    ,'wind_speed')
        s = nf90_put_att(nc,idW   , 'units'        ,'m/s')
        s = nf90_put_att(nc,idW   , 'missing_value',1.e+20)

        s = nf90_put_att(nc,idH   , 'long_name'    ,'Sea Surface Height')
        s = nf90_put_att(nc,idH   , 'units'        ,'m')
        s = nf90_put_att(nc,idH   , 'missing_value',1.e+20)

        s = nf90_put_att(nc,idE   , 'long_name'    ,'concentration/dilution water flux')
        s = nf90_put_att(nc,idE   , 'units'        ,'Kg/m2/s')
        s = nf90_put_att(nc,idE   , 'missing_value',1.e+20)


        s =nf90_enddef(nc)

        counter=0

        s = nf90_put_var(nc, idlamt,       X  ); call handle_err1(s,counter)
        s = nf90_put_var(nc, idphit,       Y  ); call handle_err1(s,counter)
        s = nf90_put_var(nc, idgdept,      Z  ); call handle_err1(s,counter)

        s = nf90_put_var(nc, idvartime, TheTime  );
        call handle_err1(s,counter)

        s= nf90_close(nc)
        PREPARE_BOX_FILE_T = .true.
    END FUNCTION PREPARE_BOX_FILE_T

