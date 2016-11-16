 program ForcingGenerator

    USE modulo16
    USE netcdf
    USE calendar


    IMPLICIT NONE

    LOGICAL, ALLOCATABLE, DIMENSION(:,:,:) :: Btmask16, Btmask8
    LOGICAL, ALLOCATABLE, DIMENSION(:,:,:) :: Bumask16, Bumask8
    LOGICAL, ALLOCATABLE, DIMENSION(:,:,:) :: Bvmask16, Bvmask8

    integer        numwri_U,  numwri_V,  numwri_W,  numwri_T  ! units of text files
    character*80   filnam_U,  filnam_V,  filnam_W,  filnam_T  ! names of text files
    character*1024 filnamU,   filnamV,   filnamW  , filnamT   ! strings inside text files (.nc files)
    character*1024 fileOUT_U, fileOUT_V, fileOUT_W, fileOUT_T ! outputs .nc

    integer mycount, ntimes

    character*100 dim_name
    integer ji,jj,jk,jt
    integer ji16,jj16,ji8,jj8      ! dimensioni mesh interpolazione
    integer im, jm, km             ! dimensioni mesh finale

    ! --------------------- usate in lettura mask_16 e mask_8
    integer xid,yid,depid, idtim
    integer idlamt,idlamu,idlamv,idlamf
    integer idphit,idphiu,idphiv,idphif
    integer ide1t,ide2t,ide1v,ide2u
    integer ide3t,ide3u,ide3v,ide3w
    integer idtmask,idumask,idvmask
    integer idgdept,idgdepw
    integer c4(4),s4(4), c3(3), s3(3),s2(2),c2(2), time

    integer vert_dyn,istr,opa_oper
    integer jpi_test,jpj_test,jpk_test
        
    REAL(4) UNDEF  ! it was Read by T16, now imposed
    logical B

!    ------ SEZIONE MAIN ------------------------------
     character*1024 ORIG_DIR, OUT_DIR, EVAPO_DIR
     character DATA17*17, UnitsTime*50, charMONTH*2
     character*1024 MASKFILE_OUT, MASKFILE_16, MASKFILE_8_INT
     character*1024 evapofile
     integer theCUT
     integer year, month, day
     double precision sec, TheTime
!   ------------ SEZIONE PS ----------------------
    logical PS
    character*1024 PS16FILE_U,PS16FILE_V
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: PS16_U, PS16_V

!   ------------  SEZIONE CUTandNUDG -------------------------

    character*1024 NUDGFILE_U,NUDGFILE_V,NUDGFILE_W,NUDGFILE_T
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: NUDG16_U, NUDG16_V, NUDG16_W, M
    REAL, ALLOCATABLE, DIMENSION(:,:)   :: NUDG16_T, X,Y,M2
    REAL, ALLOCATABLE, DIMENSION(:)     :: Z

    NAMELIST/INPUTS/ MASKFILE_16,MASKFILE_8_INT,MASKFILE_OUT,ORIG_DIR,OUT_DIR,&
                     PS,PS16FILE_U,PS16FILE_V, &
                     NUDGFILE_U,NUDGFILE_V,NUDGFILE_W,NUDGFILE_T, &
                     EVAPO_DIR

    OPEN(201, file='namelist')
    REWIND(201)
    READ(201,INPUTS)
    CLOSE(201)

    UNDEF = 1.0e+20

    UnitsTime = 'seconds since 1582-10-15T00:00:00Z'

    vert_dyn    = 1
    opa_oper    = 0 ! Only for V0
    istr        = 1


    mycount =0

    numwri_T = 292 ; filnam_T = 'nomefile_T'
    numwri_U = 293 ; filnam_U = 'nomefile_U'
    numwri_V = 294 ; filnam_V = 'nomefile_V'
    numwri_W = 295 ; filnam_W = 'nomefile_W'


    call COUNTLINE(filnam_U , ntimes)
    write(*,*) ntimes , ' time frames will be interpolated'


    B = read_mask_16()
    B = read_mask_out();
    B = logicaltmask()
    theCUT = jpi16 - im + 1

    write(*,*) 'theCUT = ', theCUT



    !------------------------------------------------
    write(*,*)'Input domain dimensions'
    write(*,*)'jpi',jpi16,'jpj',jpj16,'jpk',jpk16
   ! write(*,*)'Interp domain dimensions'
   ! write(*,*)'jpi',jpi8,'jpj',jpj8,'jpk',jpk8
    write(*,*)'Output domain dimensions'
    write(*,*)'jpi',im,'jpj',jm,'jpk',km

    !--------   Reading PS 16 -----------------------------------
    if (PS) then
        allocate(PS16_U(jpi16, jpj16, jpk16))
        allocate(PS16_V(jpi16, jpj16, jpk16))
        call readNetCDF_3dvar(PS16FILE_U,'vozocrtx',jpj16,jpi16,jpk16, PS16_U  )
        call readNetCDF_3dvar(PS16FILE_V,'vomecrty',jpj16,jpi16,jpk16, PS16_V  )
    end if

    !--------   Reading Nudg 16 -----------------------------------
    allocate(M      (im, jm, km), M2(im,jm))
    allocate(NUDG16_U(im, jm, km))
    allocate(NUDG16_V(im, jm, km))
    allocate(NUDG16_W(im, jm, km))
    allocate(NUDG16_T(im, jm    ))

    call readNetCDF_3dvar(NUDGFILE_U,'vozocrtx',im,jm,km, NUDG16_U); write(*,*) 'nudgU read'
    call readNetCDF_3dvar(NUDGFILE_V,'vomecrty',im,jm,km, NUDG16_V); write(*,*) 'nudgV read'
    call readNetCDF_3dvar(NUDGFILE_W,'vovecrtz',im,jm,km, NUDG16_W); write(*,*) 'nudgW read'
    call readNetCDF_2dvar(NUDGFILE_T,'soshfldo',im,jm,    NUDG16_T); write(*,*) 'nudgT read'


!  *************************************************************************

    open(numwri_T,file=filnam_T,form='formatted')
    open(numwri_U,file=filnam_U,form='formatted')
    open(numwri_V,file=filnam_V,form='formatted')
    open(numwri_W,file=filnam_W,form='formatted')



    do jt=1,ntimes

        write(*,*) 'Time slice n:',jt
        write(*,*) '++++++++++++++++++++++++++++++++++++++'

        ! if(NetCDF_data .EQ. 0) B=DIMG_READ();

        !------------------------------------------------
        !       Read data from MED16 (T S U V) in Netcdf Format

        read(numwri_U,'(a)') filnamU ;
        read(numwri_V,'(a)') filnamV ;
        read(numwri_W,'(a)') filnamW ;
        read(numwri_T,'(a)') filnamT ;

        DATA17 = filnamU(1:8)//'-12:00:00'
        call read_date_string(DATA17, year, month, day, sec)
        CALL ymds2ju(year,month,day,sec,TheTime)
        TheTime=TheTime*86400.0
        charMONTH = DATA17(5:6)

        fileOUT_U = trim(OUT_DIR)//'U'//DATA17//'.nc'
        fileOUT_V = trim(OUT_DIR)//'V'//DATA17//'.nc'
        fileOUT_W = trim(OUT_DIR)//'W'//DATA17//'.nc'
        fileOUT_T = trim(OUT_DIR)//'T'//DATA17//'.nc'

        filnamU = trim(ORIG_DIR)//filnamU; write(*,'(A)') 'filnamU '//trim(filnamU)
        filnamV = trim(ORIG_DIR)//filnamV;
        filnamW = trim(ORIG_DIR)//filnamW;
        filnamT = trim(ORIG_DIR)//filnamT;

        evapofile= trim(EVAPO_DIR)//'sowaflcd_yyyy'//charMONTH//'15.nc'

        sn16  =0.0
        tn16  =0.0
        un16  =0.0
        vn16  =0.0
        wn16  =0.0
        qsr16 =0.0

        B = READ_U16(); if (PS) un16 = un16*PS16_U;
        B = READ_V16(); if (PS) vn16 = vn16*PS16_V;
        B = READ_T16();

        if(vert_dyn .EQ. 1)  B = READ_W16();

        !-----------------------------------------------

        if(istr.eq.1) B = WIND_STRESS()

        !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

        B = COMPUTE_W()

!========================================
!-Writing NETCDF output






      B = PREPARE_BOX_FILE_U(trim(fileOUT_U))
      B = PREPARE_BOX_FILE_V(trim(fileOUT_V))
      B = PREPARE_BOX_FILE_W(trim(fileOUT_W))
      B = PREPARE_BOX_FILE_T(trim(fileOUT_T))




      M =  un16(theCUT:jpi16,:,:)* NUDG16_U; call MODIFY_NC_3D(trim(fileOUT_U),'vozocrtx',im,jm,km,M )
      M =  vn16(theCUT:jpi16,:,:)* NUDG16_V; call MODIFY_NC_3D(trim(fileOUT_V),'vomecrty',im,jm,km,M )

      M =  wn16(theCUT:jpi16,:,:)* NUDG16_W; call MODIFY_NC_3D(trim(fileOUT_W),'vovecrtz',im,jm,km,M )
      M = avt16(theCUT:jpi16,:,:)* NUDG16_W; call MODIFY_NC_3D(trim(fileOUT_W),'votkeavt',im,jm,km,M )

      M =  tn16(theCUT:jpi16,:,:)          ; call MODIFY_NC_3D(trim(fileOUT_T),'votemper',im,jm,km,M )
      M =  sn16(theCUT:jpi16,:,:)          ; call MODIFY_NC_3D(trim(fileOUT_T),'vosaline',im,jm,km,M )
      M2= qsr16(theCUT:jpi16,:  )* NUDG16_T; call MODIFY_NC_2D(trim(fileOUT_T),'soshfldo',im,jm,   M2)
      M2= wsp16(theCUT:jpi16,:  )* NUDG16_T; call MODIFY_NC_2D(trim(fileOUT_T),'sowindsp',im,jm,   M2)
      M2= ssh16(theCUT:jpi16,:  )          ; call MODIFY_NC_2D(trim(fileOUT_T),'sossheig',im,jm,   M2)
      M2= flp16(theCUT:jpi16,:  )          ; call MODIFY_NC_2D(trim(fileOUT_T),'sowaflcd',im,jm,   M2)

!++++++++++++++++++++++++++++++++++++++++++++++++
        enddo !t loop


       close(numwri_T)
       close(numwri_U)
       close(numwri_V)
       close(numwri_W)



        CONTAINS
! **************************************************************************
! **************************************************************************
! **************************************************************************

    LOGICAL FUNCTION READ_MASK_16()

        USE MODULO16
        IMPLICIT NONE

        integer ncid, s, ccc
        integer c2(2), s2(2)
        integer cut

        cut         = 0
        mycount     = 0
        s = nf90_open(trim(MASKFILE_16), NF90_NOWRITE, ncid) ; call handle_err1(s,mycount)

        !------------------------------------------------
        !   Get dimensios of the model (MED16) data
        !

        s = nf90_inq_dimid (ncid, 'x', xid  );  call handle_err1(s,mycount)
        s = nf90_inq_dimid (ncid, 'y', yid   ); call handle_err1(s,mycount)
        s = nf90_inq_dimid (ncid, 'z', depid ); call handle_err1(s,mycount)
        s = nf90_inq_dimid (ncid,'t',idtim); call handle_err1(s,mycount)

        s = nf90_Inquire_Dimension(ncid, xid  , dim_name, jpi16); call handle_err1(s,mycount)
        s = nf90_Inquire_Dimension(ncid, yid  , dim_name, jpj16); call handle_err1(s,mycount)
        s = nf90_Inquire_Dimension(ncid, depid, dim_name, jpk16); call handle_err1(s,mycount)
        s = nf90_Inquire_Dimension(ncid, idtim, dim_name, time ); call handle_err1(s,mycount)

        !------------------------------------------------
        !   Allocate arrays
        !

        call alloc_tot_16()

        allocate(Btmask16(jpj16,jpi16,jpk16))
        allocate(Bumask16(jpj16,jpi16,jpk16))
        allocate(Bvmask16(jpj16,jpi16,jpk16))


        !------------------------------------------------
        !   Read space coordinate values (MED16 model)
        !
        s = nf90_inq_varid (ncid, 'glamt', idlamt) ; call handle_err1(s,mycount)
        s = nf90_inq_varid (ncid, 'glamu', idlamu) ; call handle_err1(s,mycount)
        s = nf90_inq_varid (ncid, 'glamv', idlamv) ; call handle_err1(s,mycount)
        s = nf90_inq_varid (ncid, 'glamf', idlamf) ; call handle_err1(s,mycount)
        s = nf90_inq_varid (ncid, 'gphit', idphit) ; call handle_err1(s,mycount)
        s = nf90_inq_varid (ncid, 'gphiu', idphiu) ; call handle_err1(s,mycount)
        s = nf90_inq_varid (ncid, 'gphiv', idphiv) ; call handle_err1(s,mycount)
        s = nf90_inq_varid (ncid, 'gphif', idphif) ; call handle_err1(s,mycount)
        s3=(/1, 1, time/)
        c3=(/jpi16, jpj16, 1/)

        s = nf90_get_var(ncid,idlamt,glamt16,s3,c3);call handle_err1(s,mycount)
        s = nf90_get_var(ncid,idlamu,glamu16,s3,c3);call handle_err1(s,mycount)
        s = nf90_get_var(ncid,idlamv,glamv16,s3,c3);call handle_err1(s,mycount)
        s = nf90_get_var(ncid,idlamf,glamf16,s3,c3);call handle_err1(s,mycount)
        s = nf90_get_var(ncid,idphit,gphit16,s3,c3);call handle_err1(s,mycount)
        s = nf90_get_var(ncid,idphiu,gphiu16,s3,c3);call handle_err1(s,mycount)
        s = nf90_get_var(ncid,idphiv,gphiv16,s3,c3);call handle_err1(s,mycount)
        s = nf90_get_var(ncid,idphif,gphif16,s3,c3);call handle_err1(s,mycount)

        !   Read scale factors (MED16 model)
        s = nf90_inq_varid (ncid, 'e1t', ide1t); call handle_err1(s,mycount)
        s = nf90_inq_varid (ncid, 'e2t', ide2t); call handle_err1(s,mycount)
        s = nf90_inq_varid (ncid, 'e3t', ide3t); call handle_err1(s,mycount)
        s = nf90_inq_varid (ncid, 'e3u', ide3u); call handle_err1(s,mycount)
        s = nf90_inq_varid (ncid, 'e3v', ide3v); call handle_err1(s,mycount)
        s = nf90_inq_varid (ncid, 'e3w', ide3w); call handle_err1(s,mycount)
        s = nf90_inq_varid (ncid, 'e1v', ide1v); call handle_err1(s,mycount)
        s = nf90_inq_varid (ncid, 'e2u', ide2u); call handle_err1(s,mycount)
        s = nf90_inq_varid (ncid, 'tmask', idtmask) ; call handle_err1(s,mycount)
        s = nf90_inq_varid (ncid, 'umask', idumask) ; call handle_err1(s,mycount)
        s = nf90_inq_varid (ncid, 'vmask', idvmask) ; call handle_err1(s,mycount)

        s3=(/1, 1, time/)
        c3=(/jpi16, jpj16, 1/)
        s = nf90_get_var(ncid,ide1t,e1t16,s3,c3); call handle_err1(s,mycount)
        s = nf90_get_var(ncid,ide2t,e2t16,s3,c3); call handle_err1(s,mycount)
        s = nf90_get_var(ncid,ide1v,e1v16,s3,c3); call handle_err1(s,mycount)
        s = nf90_get_var(ncid,ide2u,e2u16,s3,c3); call handle_err1(s,mycount)

        s4=(/1, 1, 1, time/)
        c4=(/jpi16, jpj16, jpk16, 1/)
        s = nf90_get_var(ncid,ide3t,e3t16,s4,c4); call handle_err1(s,mycount)
        s = nf90_get_var(ncid,ide3u,e3u16,s4,c4); call handle_err1(s,mycount)
        s = nf90_get_var(ncid,ide3v,e3v16,s4,c4); call handle_err1(s,mycount)
        s = nf90_get_var(ncid,ide3w,e3w16,s4,c4); call handle_err1(s,mycount)

        s = nf90_get_var (ncid,idtmask,tmask16,s4, c4) ; call handle_err1(s,mycount)
        s = nf90_get_var (ncid,idumask,umask16,s4, c4) ; call handle_err1(s,mycount)
        s = nf90_get_var (ncid,idvmask,vmask16,s4, c4) ; call handle_err1(s,mycount)



        ! Correction due  to incongruence between meshmask and
        ! model output
        if(opa_oper .EQ. 1) then
            tmask16(199-cut,92,30)=0.
            tmask16(199-cut,93,30)=0.
            tmask16(199-cut,92,31)=0.
            tmask16(199-cut,93,31)=0.
            tmask16(199-cut,92,32)=0.
            tmask16(199-cut,93,32)=0.
            tmask16(199-cut,92,33)=0.
            tmask16(199-cut,93,33)=0.
            tmask16(199-cut,92,34)=0.
            tmask16(199-cut,93,34)=0.
        endif

        if(opa_oper .EQ. 1) then
            umask16(198-cut,92,30)=0.
            umask16(199-cut,92,30)=0.
            umask16(198-cut,93,30)=0.
            umask16(199-cut,93,30)=0.
            umask16(198-cut,92,31)=0.
            umask16(199-cut,92,31)=0.
            umask16(198-cut,93,31)=0.
            umask16(199-cut,93,31)=0.
            umask16(198-cut,92,32)=0.
            umask16(199-cut,92,32)=0.
            umask16(198-cut,93,32)=0.
            umask16(199-cut,93,32)=0.
            umask16(198-cut,92,33)=0.
            umask16(199-cut,92,33)=0.
            umask16(198-cut,93,33)=0.
            umask16(199-cut,93,33)=0.
            umask16(198-cut,92,34)=0.
            umask16(199-cut,92,34)=0.
            umask16(198-cut,93,34)=0.
            umask16(199-cut,93,34)=0.
        endif

        if(opa_oper .EQ. 1) then
            vmask16(199-cut,92,30)=0.
            vmask16(199-cut,92,31)=0.
            vmask16(199-cut,92,32)=0.
            vmask16(199-cut,92,33)=0.
            vmask16(199-cut,92,34)=0.
        endif


        c2=(/jpk16,  1/)
        s2=(/1, 1 /)


        s = nf90_inq_varid (ncid, 'gdept_0', idgdept);call handle_err1(s,mycount)
        s = nf90_inq_varid (ncid, 'gdepw_0', idgdepw); call handle_err1(s,mycount)

        s = nf90_get_var (ncid,idgdept,gdept16,s2, c2) ; call handle_err1(s,mycount)
        s = nf90_get_var (ncid,idgdepw,gdepw16,s2, c2) ; call handle_err1(s,mycount)

        s= nf90_close(ncid) ; call handle_err1(s,mycount)
        !------------------------------------------------


        READ_MASK_16 = .true.
    END FUNCTION READ_MASK_16



    ! **************************************************************************
    ! **************************************************************************

! **************************************************************************
! **************************************************************************

    LOGICAL FUNCTION READ_MASK_out()

        IMPLICIT NONE

        call getDimension(trim(MASKFILE_OUT), 'x', im)
        call getDimension(trim(MASKFILE_OUT), 'y', jm)
        call getDimension(trim(MASKFILE_OUT), 'z', km)

        allocate(X(im,jm), Y(im,jm), Z(km))

        call readNetCDF_2dvar(trim(MASKFILE_OUT), 'nav_lon',im,jm,X  );
        call readNetCDF_2dvar(trim(MASKFILE_OUT), 'nav_lat',im,jm,Y  );
        call readNetCDF_1dvar(trim(MASKFILE_OUT), 'nav_lev',km   ,Z  );


        READ_MASK_out = .true.
    END FUNCTION READ_MASK_out

! **************************************************************************
! **************************************************************************

    LOGICAL FUNCTION LOGICALTMASK()
        USE MODULO16
        IMPLICIT NONE

        Btmask16=.false.


        do ji=1,jpi16
            do jj=1,jpj16
                do jk=1,jpk16
                    if (tmask16(ji,jj,jk).eq.1)  Btmask16(ji,jj,jk)=.true.
                    if (umask16(ji,jj,jk).eq.1)  Bumask16(ji,jj,jk)=.true.
                    if (vmask16(ji,jj,jk).eq.1)  Bvmask16(ji,jj,jk)=.true.
                enddo
            enddo
        enddo

        LOGICALTMASK=.true.

    END FUNCTION LOGICALTMASK

! **************************************************************************
! **************************************************************************


    LOGICAL FUNCTION READ_U16()
        USE MODULO16
        IMPLICIT NONE

        integer nciu, s, idun, idtx

        s = nf90_open(filnamU, NF90_NOWRITE, nciu) ; call handle_err1(s,mycount)
        s = nf90_inq_dimid (nciu, 'x', xid)        ; call handle_err1(s,mycount)
        s = nf90_inq_dimid (nciu, 'y', yid)        ; call handle_err1(s,mycount)
        s = nf90_inq_dimid (nciu, 'depthu', depid) ; call handle_err1(s,mycount)

        s = nf90_Inquire_Dimension (nciu, xid  , dim_name, jpi_test) ; call handle_err1(s,mycount)
        s = nf90_Inquire_Dimension (nciu, yid  , dim_name, jpj_test) ; call handle_err1(s,mycount)
        s = nf90_Inquire_Dimension (nciu, depid, dim_name, jpk_test) ; call handle_err1(s,mycount)


        if((jpi_test-jpi16) .NE. 0 ) then
            write(*,*)'Merge_Error: dimension mismatch:meshmask,Vdata'
            write(*,*)'jpiMeshmask',jpi16, 'jpi_Udata',jpi_test
            STOP
        endif
        if((jpj_test-jpj16) .NE. 0 ) then
            write(*,*)'Merge_Error: dimension mismatch:meshmask,Udata'
            write(*,*)'jpjMeshmask',jpj16, 'jpj_Udata',jpj_test
            STOP
        endif
        if((jpk_test-jpk16) .NE. 0 ) then
            write(*,*)'Merge_Error: dimension mismatch:meshmask,Udata'
            write(*,*)'jpkMeshmask',jpk16, 'jpk_Udata',jpk_test
            STOP
        endif


        s = nf90_inq_varid (nciu, 'vozocrtx' , idun) ; call handle_err1(s,mycount)
        c4=(/jpi16, jpj16, jpk16, 1/)
        s4=(/1, 1, 1, time/)
        s = nf90_get_var (nciu, idun, un16, s4, c4) ; call handle_err1(s,mycount)


        if(istr .EQ. 1) then
            c3=(/jpi16, jpj16, 1/)
            s3=(/1, 1, time/)
            s = nf90_inq_varid (nciu, 'sozotaux' , idtx) ; call handle_err1(s,mycount)
            s = nf90_get_var (nciu, idtx, taux16, s3, c3); call handle_err1(s,mycount)
        endif

        s= nf90_close(nciu) ; call handle_err1(s,mycount)
        READ_U16 = .true.
    END FUNCTION READ_U16


! **************************************************************************
! **************************************************************************

    LOGICAL FUNCTION READ_V16()
        USE MODULO16
        IMPLICIT NONE

        integer nciv, s, idvn, idty
        s = nf90_open(filnamV, NF90_NOWRITE, nciv) ; call handle_err1(s,mycount)
        !       V data
        s = nf90_inq_dimid (nciv, 'x'     ,   xid) ; call handle_err1(s,mycount)
        s = nf90_inq_dimid (nciv, 'y'     ,   yid) ; call handle_err1(s,mycount)
        s = nf90_inq_dimid (nciv, 'depthv', depid) ; call handle_err1(s,mycount)
        s = nf90_Inquire_Dimension (nciv, xid  , dim_name, jpi_test); call handle_err1(s,mycount)
        s = nf90_Inquire_Dimension (nciv, yid  , dim_name, jpj_test); call handle_err1(s,mycount)
        s = nf90_Inquire_Dimension (nciv, depid, dim_name, jpk_test); call handle_err1(s,mycount)

        if((jpi_test-jpi16) .NE. 0 ) then
            write(*,*)'Merge_Error: dimension mismatch:meshmask,Vdata'
            write(*,*)'jpiMeshmask',jpi16, 'jpi_Vdata',jpi_test
            STOP
        endif
        if((jpj_test-jpj16) .NE. 0 ) then
            write(*,*)'Merge_Error: dimension mismatch:meshmask,Vdata'
            write(*,*)'jpjMeshmask',jpj16, 'jpj_Vdata',jpj_test
            STOP
        endif
        if((jpk_test-jpk16) .NE. 0 ) then
            write(*,*)'Merge_Error: dimension mismatch:meshmask,Vdata'
            write(*,*)'jpkMeshmask',jpk16, 'jpk_Vdata',jpk_test
            STOP
        endif

        s = nf90_inq_varid (nciv, 'vomecrty' , idvn) ; call handle_err1(s,mycount)
        c4=(/jpi16, jpj16, jpk16, 1/)
        s4=(/1, 1, 1, time/)
        s = nf90_get_var (nciv, idvn, vn16, s4, c4) ; call handle_err1(s,mycount)

        if(istr .EQ. 1) then

            s = nf90_inq_varid (nciv, 'sometauy' , idty) ; call handle_err1(s,mycount)

            c3=(/jpi16, jpj16, 1/)
            s3=(/1, 1, time/)
            s = nf90_get_var (nciv, idty, tauy16, s3, c3) ; call handle_err1(s,mycount)
        endif

        s= nf90_close(nciv) ; call handle_err1(s,mycount)
        READ_V16 = .true.
    END FUNCTION READ_V16


! **************************************************************************
! **************************************************************************

    LOGICAL FUNCTION READ_W16()
        USE MODULO16
        IMPLICIT NONE

        integer nciw, s, idwn, idavt
        s = nf90_open(filnamW, NF90_NOWRITE, nciw) ; call handle_err1(s,mycount)

        !------------------------------------------------
        !       Test on input data
        !       Get dimensios of the model (MED16) data
        !


        !       W data

        s = nf90_inq_dimid (nciw, 'x', xid)
        s = nf90_inq_dimid (nciw, 'y', yid)
        s = nf90_inq_dimid (nciw, 'depthw', depid)

        s = nf90_Inquire_Dimension (nciw, xid, dim_name, jpi_test)
        s = nf90_Inquire_Dimension (nciw, yid, dim_name, jpj_test)
        s = nf90_Inquire_Dimension (nciw, depid, dim_name, jpk_test)

        if((jpi_test-jpi16) .NE. 0 ) then
            write(*,*)'Merge_Error: dimension mismatch:meshmask,Wdata'
            write(*,*)'jpiMeshmask',jpi16, 'jpi_Wdata',jpi_test
            STOP
        endif
        if((jpj_test-jpj16) .NE. 0 ) then
            write(*,*)'Merge_Error: dimension mismatch:meshmask,Wdata'
            write(*,*)'jpjMeshmask',jpj16, 'jpj_Wdata',jpj_test
            STOP
        endif
        if((jpk_test-jpk16) .NE. 0 ) then
            write(*,*)'Merge_Error: dimension mismatch:meshmask,Wdata'
            write(*,*)'jpkMeshmask',jpk16, 'jpk_Wdata',jpk_test
            STOP
        endif

        !        s = nf90_inq_varid (nciw, 'vovecrtz' , idwn)
        s = nf90_inq_varid (nciw, 'votkeavt' , idavt) ; call handle_err1(s,mycount)

        c4=(/jpi16, jpj16, jpk16, 1/)
        s4=(/1, 1, 1, time/)

        s = nf90_get_var (nciw, idavt, avt16, s4, c4); call handle_err1(s,mycount)

        s= nf90_close(nciw) ; call handle_err1(s,mycount)

        READ_W16 = .true.
    END FUNCTION READ_W16


! **************************************************************************
! **************************************************************************


    LOGICAL FUNCTION READ_T16()
        USE MODULO16
        IMPLICIT NONE

        integer ncit, s
        integer idsow,idsoh,idsof,idssh,idsal,idtem,idflp
        s = nf90_open(filnamT, NF90_NOWRITE, ncit) ; call handle_err1(s,mycount)

        s = nf90_inq_dimid (ncit, 'x', xid) ; call handle_err1(s,mycount)
        s = nf90_inq_dimid (ncit, 'y', yid) ; call handle_err1(s,mycount)
        s = nf90_inq_dimid (ncit, 'deptht', depid) ; call handle_err1(s,mycount)

        s = nf90_Inquire_Dimension (ncit, xid, dim_name, jpi_test) ; call handle_err1(s,mycount)
        s = nf90_Inquire_Dimension (ncit, yid, dim_name, jpj_test) ; call handle_err1(s,mycount)
        s = nf90_Inquire_Dimension (ncit, depid, dim_name, jpk_test) ; call handle_err1(s,mycount)


        if((jpi_test-jpi16) .NE. 0 ) then
            write(*,*)'Merge_Error: dimension mismatch:meshmask,Tdata'
            write(*,*)'jpiMeshmask',jpi16, 'jpi_Tdata',jpi_test
            STOP
        endif
        if((jpj_test-jpj16) .NE. 0 ) then
            write(*,*)'Merge_Error: dimension mismatch:meshmask,Tdata'
            write(*,*)'jpjMeshmask',jpj16, 'jpj_Tdata',jpj_test
            STOP
        endif
        if((jpk_test-jpk16) .NE. 0 ) then
            write(*,*)'Merge_Error: dimension mismatch:meshmask,Tdata'
            write(*,*)'jpkMeshmask',jpk16, 'jpk_Tdata',jpk_test
            STOP
        endif

        s = nf90_inq_varid (ncit, 'vosaline', idsal); call handle_err1(s,mycount)
        s = nf90_inq_varid (ncit, 'votemper', idtem); call handle_err1(s,mycount)
        s = nf90_inq_varid (ncit, 'soshfldo' ,idsof); call handle_err1(s,mycount)
        s = nf90_inq_varid (ncit, 'sossheig' ,idssh); call handle_err1(s,mycount)
        s = nf90_inq_varid (ncit, 'sowaflup' ,idflp); call handle_err1(s,mycount)
        !       s = nf90_inq_varid (ncit, 'sowaflup' , idsow)
        !       s = nf90_inq_varid (ncit, 'sohefldo' , idsoh)


        c4=(/jpi16, jpj16, jpk16, 1/)
        s4=(/1, 1, 1, time/)
        s = nf90_get_var (ncit, idsal, sn16, s4,c4) ; call handle_err1(s,mycount)
        s = nf90_get_var (ncit, idtem, tn16, s4,c4) ; call handle_err1(s,mycount)

        !s = nf90_get_att (ncit, idsal,'missing_value',UNDEF) ; call handle_err1(s,mycount)
        c3=(/jpi16, jpj16, 1/)
        s3=(/1, 1, time/)
        s = nf90_get_var (ncit, idsof, qsr16, s3,c3) ; call handle_err1(s,mycount)
        s = nf90_get_var (ncit, idssh, ssh16, s3,c3) ; call handle_err1(s,mycount)
        s = nf90_get_var (ncit, idflp, flp16, s3,c3) ; call handle_err1(s,mycount)

        s= nf90_close(ncit) ; call handle_err1(s,mycount)
        READ_T16 = .true.
    END FUNCTION READ_T16


    ! **************************************************************************
    ! **************************************************************************


    LOGICAL FUNCTION WIND_STRESS()
        USE MODULO16
        IMPLICIT NONE

        ! Velocity definition using wind stress
        !      very rude treatement of wind stress extrapolation
        !      i.e. where (tmask16 .EQ. 1) .AND. (umask16 .EQ.0)

        taux16(1,:)=0
        tauy16(1,:)=0
        wsp16(1,:)=0

        taux16(:,1)=0
        tauy16(:,1)=0
        wsp16(:,1)=0

        taux16(jpi16,:)=0
        tauy16(jpi16,:)=0
        wsp16(jpi16,:)=0

        taux16(:,jpj16)=0
        tauy16(:,jpj16)=0
        wsp16(:,jpj16)=0
        aux16=taux16
        auy16=tauy16

        do jj=2,jpj16
            do ji=2,jpi16

                if (Btmask16(ji,jj,1).AND.(.not.Bumask16(ji,jj,1))) then
                    aux16_01(ji,jj) = min(taux16(ji+1,jj),taux16(ji+1,jj+1))
                    auy16_01(ji,jj) = min(tauy16(ji+1,jj),tauy16(ji+1,jj+1))

                    aux16_02(ji,jj) = min(aux16_01(ji,jj),taux16(ji,jj+1))
                    auy16_02(ji,jj) = min(auy16_01(ji,jj),tauy16(ji,jj+1))

                    aux16_03(ji,jj) = min(aux16_02(ji,jj),taux16(ji-1,jj+1))
                    auy16_03(ji,jj) = min(auy16_02(ji,jj),tauy16(ji-1,jj+1))

                    aux16_04(ji,jj) = min(aux16_03(ji,jj),taux16(ji-1,jj))
                    auy16_04(ji,jj) = min(auy16_03(ji,jj),tauy16(ji-1,jj))

                    aux16_05(ji,jj) = min(aux16_04(ji,jj),taux16(ji-1,jj-1))
                    auy16_05(ji,jj) = min(auy16_04(ji,jj),tauy16(ji-1,jj-1))

                    aux16_06(ji,jj) = min(aux16_05(ji,jj),taux16(ji,jj-1))
                    auy16_06(ji,jj) = min(auy16_05(ji,jj),tauy16(ji,jj-1))

                    aux16(ji,jj) = min(aux16_06(ji,jj),taux16(ji+1,jj-1))
                    auy16(ji,jj) = min(auy16_06(ji,jj),tauy16(ji+1,jj-1))
                endif

                if(Btmask16(ji,jj,1) .AND. (.not.Bvmask16(ji,jj,1))) then
                    aux16_01(ji,jj) = min(taux16(ji+1,jj),taux16(ji+1,jj+1))
                    auy16_01(ji,jj) = min(tauy16(ji+1,jj),tauy16(ji+1,jj+1))

                    aux16_02(ji,jj) = min(aux16_01(ji,jj),taux16(ji,jj+1))
                    auy16_02(ji,jj) = min(auy16_01(ji,jj),tauy16(ji,jj+1))

                    aux16_03(ji,jj) = min(aux16_02(ji,jj),taux16(ji-1,jj+1))
                    auy16_03(ji,jj) = min(auy16_02(ji,jj),tauy16(ji-1,jj+1))

                    aux16_04(ji,jj) = min(aux16_03(ji,jj),taux16(ji-1,jj))
                    auy16_04(ji,jj) = min(auy16_03(ji,jj),tauy16(ji-1,jj))

                    aux16_05(ji,jj) = min(aux16_04(ji,jj),taux16(ji-1,jj-1))
                    auy16_05(ji,jj) = min(auy16_04(ji,jj),tauy16(ji-1,jj-1))

                    aux16_06(ji,jj) = min(aux16_05(ji,jj),taux16(ji,jj-1))
                    auy16_06(ji,jj) = min(auy16_05(ji,jj),tauy16(ji,jj-1))

                    aux16(ji,jj) = min(aux16_06(ji,jj),taux16(ji+1,jj-1))
                    auy16(ji,jj) = min(auy16_06(ji,jj),tauy16(ji+1,jj-1))
                endif

                if(Btmask16(ji,jj,1) .AND. (.not.Bumask16(ji,jj,1)) .AND. (.not.Bvmask16(ji,jj,1))) then
                    aux16_01(ji,jj) = min(taux16(ji+1,jj),taux16(ji+1,jj+1))
                    auy16_01(ji,jj) = min(tauy16(ji+1,jj),tauy16(ji+1,jj+1))

                    aux16_02(ji,jj) = min(aux16_01(ji,jj),taux16(ji,jj+1))
                    auy16_02(ji,jj) = min(auy16_01(ji,jj),tauy16(ji,jj+1))

                    aux16_03(ji,jj) = min(aux16_02(ji,jj),taux16(ji-1,jj+1))
                    auy16_03(ji,jj) = min(auy16_02(ji,jj),tauy16(ji-1,jj+1))

                    aux16_04(ji,jj) = min(aux16_03(ji,jj),taux16(ji-1,jj))
                    auy16_04(ji,jj) = min(auy16_03(ji,jj),tauy16(ji-1,jj))

                    aux16_05(ji,jj) = min(aux16_04(ji,jj),taux16(ji-1,jj-1))
                    auy16_05(ji,jj) = min(auy16_04(ji,jj),tauy16(ji-1,jj-1))

                    aux16_06(ji,jj) = min(aux16_05(ji,jj),taux16(ji,jj-1))
                    auy16_06(ji,jj) = min(auy16_05(ji,jj),tauy16(ji,jj-1))

                    aux16(ji,jj) = min(aux16_06(ji,jj),taux16(ji+1,jj-1))
                    auy16(ji,jj) = min(auy16_06(ji,jj),tauy16(ji+1,jj-1))
                endif
            enddo
        enddo
       
        do jj=2,jpj16
            do ji=2,jpi16
                taux16(ji,jj) = aux16(ji,jj)
                tauy16(ji,jj) = auy16(ji,jj)
            enddo
        enddo

      call PURE_WIND_SPEED(taux16,tauy16,jpj16,jpi16,wsp16);

!        C1=1/(1.3*1.5*0.001)
!        do jj=2,jpj16
!            do ji=2,jpi16
!
!                wsp16(ji,jj)=sqrt(C1*sqrt(taux16(ji,jj)*taux16(ji,jj)+ &
!                tauy16(ji,jj)*tauy16(ji,jj)))*tmask16(ji,jj,1)
!            enddo
!        enddo
       
        do jj=2,jpj16
            do ji=2,jpi16
                if(.not.Btmask16(ji,jj,1)) then
                    wsp16(ji,jj)= UNDEF
                endif
            enddo
        enddo
       

        WIND_STRESS=.true.
    END FUNCTION WIND_STRESS
   ! **************************************************************************
    ! **************************************************************************



    LOGICAL FUNCTION COMPUTE_W()
        USE MODULO16
        IMPLICIT NONE

        do jk =1,jpk16
            !C-CC Zero velocities on the stencil

            un16(    1,    :,jk) = 0
            un16(    :,    1,jk) = 0
            un16(jpi16,    :,jk) = 0
            un16(    :,jpj16,jk) = 0

            vn16(    1,    :,jk) = 0
            vn16(    :,    1,jk) = 0
            vn16(jpi16,    :,jk) = 0
            vn16(    :,jpj16,jk) = 0

            hdivn16(    1,    :,jk) = 0
            hdivn16(    :,    1,jk) = 0
            hdivn16(jpi16,    :,jk) = 0
            hdivn16(    :,jpj16,jk) = 0

            do jj=1,jpj16
                do ji=1,jpi16
                    if(.not.Bumask16(ji,jj,jk)) un16(ji,jj,jk)=0
                    if(.not.Bvmask16(ji,jj,jk)) vn16(ji,jj,jk)=0
                enddo
            enddo

            ! ----------------------------------------------------------------
            call div16(jk)
            ! ----------------------------------------------------------------
            do jj=1,jpj16
                do ji=1,jpi16
                    if(.not.Bumask16(ji,jj,jk)) un16(ji,jj,jk)=0
                    if(.not.Bvmask16(ji,jj,jk)) vn16(ji,jj,jk)=0.
                enddo
            enddo


        enddo ! End on loop on jk

        ! ----------------------------------------------------------------
            call wzv16()
        ! ----------------------------------------------------------------
        do jk=1,jpk16
            do jj=1,jpj16
                do ji=1,jpi16
                    if(.not.Btmask16(ji,jj,jk)) wn16(ji,jj,jk)=0.
                enddo
            enddo
        enddo

        ! ----------------------------------------------------------------
        write(*,*) 'comp_ssh16'
        call comp_ssh16
        ! ----------------------------------------------------------------
        do jk=1,jpk16
            call div16(jk)
        enddo

        call wzv16()

        do jk=1,jpk16
            do jj=1,jpj16
                do ji=1,jpi16
                    if(.not.Bumask16(ji,jj,jk)) un16(ji,jj,jk)=UNDEF
                    if(.not.Bvmask16(ji,jj,jk)) vn16(ji,jj,jk)=UNDEF
                    if(.not.Btmask16(ji,jj,jk)) wn16(ji,jj,jk)=UNDEF
                enddo
            enddo
        enddo
        do jj=1,jpj16
            do ji=1,jpi16
                if(.not.Btmask16(ji,jj,1)) ssh16(ji,jj)=UNDEF
            enddo
        enddo


        COMPUTE_W = .true.
    END FUNCTION COMPUTE_W




! **************************************************************************
        LOGICAL FUNCTION DIMG_READ()

            IMPLICIT NONE
            integer dimgidtx, dimgidty, dimgidQ, dimgidS, dimgidT
            integer dimgidU, dimgidV, dimgidW
            character*80 filnamtx,filnamty,filnamQ,filnamS
            character*80 filnam_2D, filnam_S
            integer     numwri_2D,numwri_S
            character    clver*4, line1*80
            write(*,*)' dimg format for data input selected'

            numwri_2D = 290
            numwri_S = 291
            dimgidtx=84
            dimgidty=86
            dimgidQ=88
            dimgidS=90
            dimgidT=92
            dimgidU=94
            dimgidV=96
            dimgidW=98
            !-------
            !       2 dimensional data
            read(numwri_2D,*) filnam_2D
            write(*,*) 'filnam_2D ', filnam_2D
            if(istr .EQ. 1) then
                write(*,*)' zonal wind stress data dimg read'
                filnamtx = filnam_2D
                open(dimgidtx,file=filnamtx,form='unformatted',access='direct',recl= jpi16*jpj16*4)
                !       Record for taux is 2
                jk=2
                read(dimgidtx,rec=jk) ((taux16(ji,jj) ,ji=1,jpi16),jj=1,jpj16)

                close(dimgidtx)
                !-------
                write(*,*)' meridional wind stress data dimg read'
                filnamty = filnam_2D
                open(dimgidty,file=filnamtx,form='unformatted',access='direct',recl= jpi16*jpj16*4)
                !       Record for tauy is 3
                jk=3
                read(dimgidty,rec=jk) ((tauy16(ji,jj) ,ji=1,jpi16),jj=1,jpj16)
                close(dimgidty)

                ! Velocity definition using wind stress
                !      very rude treatement of wind stress extrapolation
                !      where (tmask16 .EQ. 1) .AND. (umask16 .EQ.0)

                taux16(1,:)=0
                tauy16(1,:)=0
                wsp16(1,:)=0

                taux16(:,1)=0
                tauy16(:,1)=0
                wsp16(:,1)=0

                taux16(jpi16,:)=0
                tauy16(jpi16,:)=0
                wsp16(jpi16,:)=0

                taux16(:,jpj16)=0
                tauy16(:,jpj16)=0
                wsp16(:,jpj16)=0
                aux16=taux16
                auy16=tauy16

                do jj=2,jpj16
                    do ji=2,jpi16

                        if((tmask16(ji,jj,1) .EQ. 1) .AND. (umask16(ji,jj,1) .EQ.0)) then
                            aux16_01(ji,jj) = min(taux16(ji+1,jj),taux16(ji+1,jj+1))
                            auy16_01(ji,jj) = min(tauy16(ji+1,jj),tauy16(ji+1,jj+1))

                            aux16_02(ji,jj) = min(aux16_01(ji,jj),taux16(ji,jj+1))
                            auy16_02(ji,jj) = min(auy16_01(ji,jj),tauy16(ji,jj+1))

                            aux16_03(ji,jj) = min(aux16_02(ji,jj),taux16(ji-1,jj+1))
                            auy16_03(ji,jj) = min(auy16_02(ji,jj),tauy16(ji-1,jj+1))

                            aux16_04(ji,jj) = min(aux16_03(ji,jj),taux16(ji-1,jj))
                            auy16_04(ji,jj) = min(auy16_03(ji,jj),tauy16(ji-1,jj))

                            aux16_05(ji,jj) = min(aux16_04(ji,jj),taux16(ji-1,jj-1))
                            auy16_05(ji,jj) = min(auy16_04(ji,jj),tauy16(ji-1,jj-1))

                            aux16_06(ji,jj) = min(aux16_05(ji,jj),taux16(ji,jj-1))
                            auy16_06(ji,jj) = min(auy16_05(ji,jj),tauy16(ji,jj-1))

                            aux16(ji,jj) = min(aux16_06(ji,jj),taux16(ji+1,jj-1))
                            auy16(ji,jj) = min(auy16_06(ji,jj),tauy16(ji+1,jj-1))
                        endif

                        if((tmask16(ji,jj,1) .EQ. 1) .AND. (vmask16(ji,jj,1) .EQ.0)) then
                            aux16_01(ji,jj) = min(taux16(ji+1,jj),taux16(ji+1,jj+1))
                            auy16_01(ji,jj) = min(tauy16(ji+1,jj),tauy16(ji+1,jj+1))

                            aux16_02(ji,jj) = min(aux16_01(ji,jj),taux16(ji,jj+1))
                            auy16_02(ji,jj) = min(auy16_01(ji,jj),tauy16(ji,jj+1))

                            aux16_03(ji,jj) = min(aux16_02(ji,jj),taux16(ji-1,jj+1))
                            auy16_03(ji,jj) = min(auy16_02(ji,jj),tauy16(ji-1,jj+1))

                            aux16_04(ji,jj) = min(aux16_03(ji,jj),taux16(ji-1,jj))
                            auy16_04(ji,jj) = min(auy16_03(ji,jj),tauy16(ji-1,jj))

                            aux16_05(ji,jj) = min(aux16_04(ji,jj),taux16(ji-1,jj-1))
                            auy16_05(ji,jj) = min(auy16_04(ji,jj),tauy16(ji-1,jj-1))

                            aux16_06(ji,jj) = min(aux16_05(ji,jj),taux16(ji,jj-1))
                            auy16_06(ji,jj) = min(auy16_05(ji,jj),tauy16(ji,jj-1))

                            aux16(ji,jj) = min(aux16_06(ji,jj),taux16(ji+1,jj-1))
                            auy16(ji,jj) = min(auy16_06(ji,jj),tauy16(ji+1,jj-1))
                        endif

                        if((tmask16(ji,jj,1) .EQ. 1) .AND. (umask16(ji,jj,1) .EQ.0) &
                        .AND. (vmask16(ji,jj,1) .EQ.0) ) then
                            aux16_01(ji,jj) = min(taux16(ji+1,jj),taux16(ji+1,jj+1))
                            auy16_01(ji,jj) = min(tauy16(ji+1,jj),tauy16(ji+1,jj+1))

                            aux16_02(ji,jj) = min(aux16_01(ji,jj),taux16(ji,jj+1))
                            auy16_02(ji,jj) = min(auy16_01(ji,jj),tauy16(ji,jj+1))

                            aux16_03(ji,jj) = min(aux16_02(ji,jj),taux16(ji-1,jj+1))
                            auy16_03(ji,jj) = min(auy16_02(ji,jj),tauy16(ji-1,jj+1))

                            aux16_04(ji,jj) = min(aux16_03(ji,jj),taux16(ji-1,jj))
                            auy16_04(ji,jj) = min(auy16_03(ji,jj),tauy16(ji-1,jj))

                            aux16_05(ji,jj) = min(aux16_04(ji,jj),taux16(ji-1,jj-1))
                            auy16_05(ji,jj) = min(auy16_04(ji,jj),tauy16(ji-1,jj-1))

                            aux16_06(ji,jj) = min(aux16_05(ji,jj),taux16(ji,jj-1))
                            auy16_06(ji,jj) = min(auy16_05(ji,jj),tauy16(ji,jj-1))

                            aux16(ji,jj) = min(aux16_06(ji,jj),taux16(ji+1,jj-1))
                            auy16(ji,jj) = min(auy16_06(ji,jj),tauy16(ji+1,jj-1))
                        endif
                    enddo
                enddo

                do jj=2,jpj16
                    do ji=2,jpi16
                        taux16(ji,jj) = aux16(ji,jj)
                        tauy16(ji,jj) = auy16(ji,jj)
                    enddo
                enddo

                call PURE_WIND_SPEED(taux16,tauy16,jpj16,jpi16,wsp16);
!                C1=1/(1.3*1.5*0.001)
!                do jj=2,jpj16
!                    do ji=2,jpi16
!                        wsp16(ji,jj)=sqrt(C1*sqrt(taux16(ji,jj)*taux16(ji,jj)+ tauy16(ji,jj)*tauy16(ji,jj)))
!                    enddo
!                enddo

            endif ! if(istr .EQ. 1)
            !-------
            write(*,*)' qsr data dimg read'
            filnamQ = filnam_2D
            open(dimgidQ,file=filnamQ,form='unformatted', access='direct',recl= jpi16*jpj16*4)
            read(dimgidQ,rec=1) clver, line1
            !       Record for qsr is 8
            jk=8
            read(dimgidQ,rec=jk) ((qsr16(ji,jj) ,ji=1,jpi16),jj=1,jpj16)
            close(dimgidQ)

            !-------
            write(*,*)' Salinity data dimg read'
            read(numwri_S,*) filnam_S
            filnamS = filnam_S
            write(*,*) 'filnamS ', filnamS
            open(dimgidS,file=filnamS,form='unformatted', access='direct',recl= jpi16*jpj16*4)
            read(dimgidS,rec=1) clver, line1
            write(*,*) clver,line1
            do jk=1,jpk16
                read(dimgidS,rec=jk+1) ((sn16(ji,jj,jk) ,ji=1,jpi16),jj=1,jpj16)
            enddo
            close(dimgidS)

            !-------
            write(*,*)' Temperature data dimg read'
            read(numwri_T,*) filnam_T
            filnamT = filnam_T
            write(*,*) 'filnamT ', filnamT
            open(dimgidT,file=filnamT,form='unformatted', access='direct',recl= jpi16*jpj16*4)
            read(dimgidT,rec=1) clver, line1
            write(*,*) clver,line1
            do jk=1,jpk16
                read(dimgidT,rec=jk+1) ((tn16(ji,jj,jk),ji=1,jpi16),jj=1,jpj16)
            enddo
            close(dimgidT)

            !-------
            write(*,*)' Zonal velocity data dimg read'
            read(numwri_U,*) filnam_U
            filnamU = filnam_U
            write(*,*) 'filnamU ', filnamU
            open(dimgidU,file=filnamU,form='unformatted', access='direct',recl= jpi16*jpj16*4)
            read(dimgidU,rec=1) clver,line1
            write(*,*) clver,line1
            do jk=1,jpk16
                read(dimgidU,rec=jk+1) ((un16(ji,jj,jk) ,ji=1,jpi16),jj=1,jpj16)
            enddo
            close(dimgidU)

            !-------
            write(*,*)' meridional velocity data dimg read'
            read(numwri_V,*) filnam_V
            filnamV = filnam_V
            write(*,*) 'filnamV ', filnamV
            open(dimgidV,file=filnamV,form='unformatted', access='direct',recl= jpi16*jpj16*4)
            read(dimgidV,rec=1) clver, line1
            write(*,*) clver,line1
            do jk=1,jpk16
                read(dimgidV,rec=jk+1) ((vn16(ji,jj,jk) ,ji=1,jpi16),jj=1,jpj16)
            enddo
            close(dimgidV)

            if(vert_dyn .EQ. 1) then
                !-------
                write(*,*)'Vertical Eddy diffusivity  data dimg read'
                read(numwri_W,*) filnam_W
                filnamW = filnam_W
                write(*,*) 'filnamW ', filnamW
                open(dimgidW,file=filnamW,form='unformatted', access='direct',recl= jpi16*jpj16*4)
                read(dimgidW,rec=1) clver, line1
                write(*,*) clver,line1
                do jk=1,jpk16
                    read(dimgidW,rec=jk+1) ((avt16(ji,jj,jk) ,ji=1,jpi16),jj=1,jpj16)
                enddo
                close(dimgidW)
            endif
            ! UNDEF define not loaded by file in this case
            UNDEF = 1.e+20
            !-----------------------------------------------
            write(*,*)'end of dimg  file reading'



            DIMG_READ=.TRUE.
        END FUNCTION DIMG_READ

! **************************************************************************
! **************************************************************************
#include "ForcingBox.f90"
END PROGRAM ForcingGenerator

! **************************************************************************


SUBROUTINE PURE_WIND_SPEED(TAUx, TAUy, sizex,sizey, WSP)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: sizex, sizey
    REAL, DIMENSION(sizex,sizey), INTENT(IN ) :: TAUx, TAUy
    REAL, DIMENSION(sizex,sizey), INTENT(OUT) :: WSP
    ! local
    REAL rho, Cdrag, K
    rho    = 1.3 ! kg/m3
    Cdrag  = 1.5 * 0.001

    K     = sqrt(1/(rho*Cdrag));
    WSP = (TAUx**2 + TAUy**2)**0.25 * K


END SUBROUTINE PURE_WIND_SPEED

! *************************************************************************

        SUBROUTINE COUNTLINE(FILENAME,LINES)
        character FILENAME*(*)
        integer lines

        lines=0
        OPEN(UNIT=1,file=FILENAME,status='old')
        DO WHILE (.true.)
           read(1, *, END=21)
           lines = lines+1
        ENDDO

        21 CLOSE(1)

        END SUBROUTINE COUNTLINE




