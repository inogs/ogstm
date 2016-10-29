       SUBROUTINE EXISTVAR(fileNetCDF, varname,B)
       USE netcdf
       IMPLICIT NONE
       character, INTENT(IN) :: fileNetCDF*(*) ,varname*(*)
       LOGICAL,   INTENT(OUT)::B
       ! local
       integer stat,ncid,VARid,counter

       counter=0
       B=.false.
         stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)  
       call handle_err1(stat, counter,FileNetCDF)
         stat = nf90_inq_varid (ncid, varname, VARid)

         if(stat .ne. nf90_NoErr)  then
           !write(*,*) 'ERROR in Var = ', varname, ' file :', fileNetCDF
           B = .false.
         else
            B=.true.
        endif
        stat = nf90_close(ncid)                           
       call handle_err1(stat, counter,FileNetCDF)
       END SUBROUTINE EXISTVAR


! ! *****************************************************************************

      SUBROUTINE readnc_slice_double(fileNetCDF,varname, M)
      USE myalloc
      ! epascolo USE myalloc_mpp
      USE netcdf
      implicit none


      character,intent(in) :: fileNetCDF*(*) ,varname*(*)
      double precision,intent(inout) ::  M(jpk,jpj,jpi)
      
      double precision,allocatable,dimension(:,:,:) :: copy_in
      integer ncid, stat, VARid,i,j,k
      integer counter
      integer thecount(4), start(4)

      allocate(copy_in(jpi,jpj,jpk))
      counter = 0
      start    = (/nimpp, njmpp,  1,  1/)
      thecount = (/jpi,     jpj, jpk, 1/)

      stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)  
       call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_inq_varid (ncid, varname, VARid)      
       call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_get_var (ncid,VARid,copy_in,start, thecount)
      
      call handle_err2(stat, fileNetCDF,varname)        
      call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_close(ncid)                           
      call handle_err1(stat, counter,FileNetCDF)

      DO i=1,jpi
        DO j=1,jpj
          DO k=1,jpk
            M(k,j,i) = copy_in(i,j,k)
          ENDDO
        ENDDO
      ENDDO   
      
      deallocate(copy_in)

      END SUBROUTINE readnc_slice_double

! ********************************************************************************************

      SUBROUTINE readnc_slice_float(fileNetCDF,varname, M)
      USE myalloc
      ! epascolo USE myalloc_mpp
      USE netcdf
      implicit none


      character,intent(in) :: fileNetCDF*(*) ,varname*(*)
      double precision,intent(inout) ::  M(jpk,jpj,jpi)
      
      double precision,allocatable,dimension(:,:,:) :: copy_in
      integer ncid, stat, VARid,i,j,k
      integer counter
      integer thecount(4), start(4)

      allocate(copy_in(jpi,jpj,jpk))
      counter = 0
      start    = (/nimpp, njmpp,  1,  1/)
      thecount = (/jpi,     jpj, jpk, 1/)

      stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)  
       call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_inq_varid (ncid, varname, VARid)      
       call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_get_var (ncid,VARid,copy_in,start, thecount)
      
      call handle_err2(stat, fileNetCDF,varname)        
      call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_close(ncid)                           
      call handle_err1(stat, counter,FileNetCDF)

      DO i=1,jpi
        DO j=1,jpj
          DO k=1,jpk
            M(k,j,i) = real(copy_in(i,j,k),8)
          ENDDO
        ENDDO
      ENDDO   
      
      deallocate(copy_in)

      END SUBROUTINE readnc_slice_float

! ********************************************************************************************

      SUBROUTINE readnc_slice_int(fileNetCDF,varname, M)
      USE myalloc
      USE netcdf

      implicit none

      character,intent(in) :: fileNetCDF*(*) ,varname*(*)
      integer,intent(inout) ::  M(jpk,jpj,jpi)
      integer,allocatable,dimension(:,:,:) :: copy_in
      integer ncid, stat, VARid,i,j,k
      integer counter
      integer thecount(4), start(4)
      
      counter = 0
      allocate(copy_in(jpi,jpj,jpk))
      
      start    = (/nimpp, njmpp,  1,  1/)
      thecount = (/jpi,     jpj, jpk, 1/)

      stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)  
      call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_inq_varid (ncid, varname, VARid)      
      call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_get_var (ncid,VARid,copy_in,start, thecount)
      call handle_err2(stat, fileNetCDF,varname)        
      call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_close(ncid)                           
      call handle_err1(stat, counter,FileNetCDF)

      DO i=1,jpi
        DO j=1,jpj
          DO k=1,jpk
            M(k,j,i) = copy_in(i,j,k)
          ENDDO
        ENDDO
      ENDDO
       deallocate(copy_in)

      END SUBROUTINE readnc_slice_int

! ********************************************************************************************
! ********************************************************************************************
! ********************************************************************************************

      SUBROUTINE readnc_slice_double_2d(fileNetCDF,varname, M)
      USE myalloc
      ! epascolo USE myalloc_mpp
      USE netcdf
      implicit none

      character ,intent(in) :: fileNetCDF*(*) ,varname*(*)
      double precision,intent(inout) :: M(jpj,jpi)
      double precision,allocatable,dimension(:,:) :: copy_in
      
      integer ncid, stat, VARid,i,j,k
      integer counter
      integer thecount(3), start(3)

      allocate(copy_in(jpi,jpj))
      counter = 0
      start    = (/nimpp, njmpp,  1/)
      thecount = (/jpi,     jpj,  1/)


      stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)  
       call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_inq_varid (ncid, varname, VARid)      
       call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_get_var (ncid,VARid,copy_in,start, thecount)
      
      call handle_err2(stat, fileNetCDF,varname)        
       call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_close(ncid)                           
       call handle_err1(stat, counter,FileNetCDF)

      DO i=1,jpi
        DO j=1,jpj
            M(j,i) = copy_in(i,j)  
        ENDDO
      ENDDO  

      deallocate(copy_in)


      END SUBROUTINE readnc_slice_double_2d


! ********************************************************************************************
! ********************************************************************************************
! ********************************************************************************************

     SUBROUTINE readnc_slice_float_2d(fileNetCDF,varname, M)
     USE myalloc
      ! epascolo USE myalloc_mpp
      USE netcdf
      implicit none

      character ,intent(in) :: fileNetCDF*(*) ,varname*(*)
      double precision,intent(inout) :: M(jpj,jpi)
      real,allocatable,dimension(:,:) :: copy_in
      
      integer ncid, stat, VARid,i,j,k
      integer counter
      integer thecount(3), start(3)

      allocate(copy_in(jpi,jpj))
      counter = 0
      start    = (/nimpp, njmpp,  1/)
      thecount = (/jpi,     jpj,  1/)


      stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)  
       call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_inq_varid (ncid, varname, VARid)      
       call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_get_var (ncid,VARid,copy_in,start, thecount)
      
      call handle_err2(stat, fileNetCDF,varname)        
       call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_close(ncid)                           
       call handle_err1(stat, counter,FileNetCDF)

      DO i=1,jpi
        DO j=1,jpj
            M(j,i) = real(copy_in(i,j),8)  
        ENDDO
      ENDDO  

     deallocate(copy_in)


     END SUBROUTINE readnc_slice_float_2d


! ************************************************************
      SUBROUTINE readnc_global_double(fileNetCDF,varname, M)
      USE myalloc
      USE netcdf
      implicit none


      character,intent(in) :: fileNetCDF*(*) ,varname*(*)
      double precision,intent(inout) ::  M(jpk,jpjglo,jpiglo)
      
      double precision,allocatable,dimension(:,:,:) :: copy_in
      
      integer ncid, stat, VARid,i,j,k
      integer counter
      integer thecount(4), start(4)
      
      allocate(copy_in(jpiglo,jpjglo,jpk))
      counter = 0
      start    = (/1,       1,       1,  1/)
      thecount = (/jpiglo,  jpjglo,  jpk,1/)


      stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)  
       call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_inq_varid (ncid, varname, VARid)      
       call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_get_var (ncid,VARid,copy_in,start, thecount)
      call handle_err2(stat, fileNetCDF,varname)        
       call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_close(ncid)                           
       call handle_err1(stat, counter,FileNetCDF)

      DO i=1,jpiglo
       DO j=1,jpjglo
        DO k=1,jpk
          M(k,j,i) = copy_in(i,j,k)
        ENDDO
       ENDDO
      ENDDO  

      deallocate(copy_in)

      END SUBROUTINE readnc_global_double


! ************************************************************

      SUBROUTINE readnc_global_double_2d(fileNetCDF,varname, M)
      USE myalloc
      USE netcdf
      implicit none
      
      character,intent(in) :: fileNetCDF*(*) ,varname*(*)
      double precision,intent(inout) :: M(jpjglo,jpiglo)

      double precision,allocatable,dimension(:,:) :: copy_in
      integer ncid, stat, VARid,i,j
      integer counter
      integer thecount(4), start(4)

      allocate(copy_in(jpiglo,jpjglo))
      counter = 0
      
      start    = (/1,       1,      1, 1/)
! epascolo warning      
      thecount = (/jpiglo,  jpjglo, 1, 1/)

      stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)  
      call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_inq_varid (ncid, varname, VARid)      
      call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_get_var (ncid,VARid,copy_in,start, thecount)
      call handle_err2(stat, fileNetCDF,varname)        
      call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_close(ncid)                           
      call handle_err1(stat, counter,FileNetCDF)

      DO i=1,jpiglo
        DO j=1,jpjglo
          M(j,i) = copy_in(i,j)
        ENDDO
      ENDDO      

      deallocate(copy_in)

      END SUBROUTINE readnc_global_double_2d



        !****************************************************************************
        SUBROUTINE readmask_double_1d(fileNetCDF,varname,ARRAY)
        use netcdf
        USE modul_param
        implicit none
        character,intent(in) :: fileNetCDF*(*) ,varname*(*)
        double precision,intent(inout) :: ARRAY(jpk)

        integer ncid, stat, VARid
        integer counter
        integer thecount(4), start(4)

        counter=0
        start    = (/1,       1,       1,  1/)
        thecount = (/1,       1,      jpk, 1/)

        stat = nf90_open(fileNetCDF, nf90_nowrite, ncid) 
       call handle_err1(stat, counter,FileNetCDF)
        stat = nf90_inq_varid (ncid, varname, VARid)     
       call handle_err1(stat, counter,FileNetCDF)
        stat = nf90_get_var (ncid,VARid,ARRAY,start, thecount)
      
        call handle_err2(stat, fileNetCDF,varname)       
       call handle_err1(stat, counter,FileNetCDF)
        stat = nf90_close(ncid)                          
       call handle_err1(stat, counter,FileNetCDF)


        end SUBROUTINE readmask_double_1d


        !****************************************************************************
        SUBROUTINE readnc_double_1d(fileNetCDF,varname,im,ARRAY)
        use netcdf
        use myalloc
        implicit none

        character,intent(in) :: fileNetCDF*(*) ,varname*(*)
        integer,intent(in) :: im
        double precision,intent(inout),dimension(im) :: ARRAY
        
        integer ncid, stat, VARid
        integer counter
        
        counter=0

        stat = nf90_open(fileNetCDF, nf90_nowrite, ncid) 
       call handle_err1(stat, counter,FileNetCDF)
        stat = nf90_inq_varid (ncid, varname, VARid)     
       call handle_err1(stat, counter,FileNetCDF)
        stat = nf90_get_var (ncid,VARid,ARRAY)           
      
        call handle_err2(stat, fileNetCDF,varname)       
       call handle_err1(stat, counter,FileNetCDF)
        stat = nf90_close(ncid)                          
       call handle_err1(stat, counter,FileNetCDF)


        end SUBROUTINE readnc_double_1d

        !****************************************************************************
        SUBROUTINE readnc_int_1d(fileNetCDF,varname,im,ARRAY)
        use netcdf
        use myalloc
        implicit none

        character,intent(in) :: fileNetCDF*(*) ,varname*(*)
        double precision,intent(inout) :: ARRAY(jpk)
        integer,intent(in) :: im

        integer ncid, stat, VARid
        integer counter
        

        counter=0

       stat = nf90_open(fileNetCDF, nf90_nowrite, ncid) 
       call handle_err1(stat, counter,FileNetCDF)
       stat = nf90_inq_varid (ncid, varname, VARid)     
       call handle_err1(stat, counter,FileNetCDF)
       stat = nf90_get_var (ncid,VARid,ARRAY)
       call handle_err2(stat, fileNetCDF,varname)       
       call handle_err1(stat, counter,FileNetCDF)
       stat = nf90_close(ncid)                          
       call handle_err1(stat, counter,FileNetCDF)


        end SUBROUTINE readnc_int_1d

!        !****************************************************************************


        SUBROUTINE getDIMENSION(fileNetCDF,dimname,n)
        use netcdf
        implicit none

        character,intent(in) :: fileNetCDF*(*) ,dimname*(*)
        integer,intent(inout) :: n
       
        

        ! local

        integer DIMid,ncid,stat
        character(LEN=100) junk
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
        END SUBROUTINE getDIMENSION

       !****************************************************************************
       !****************************************************************************
       !****************************************************************************

       SUBROUTINE write_restart(fileNetCDF,VAR, julian)
       USE netcdf
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
!      istart=index(fileNetCDF,'00.')+3
!      iend  =index(fileNetCDF,'.nc')-1
!      VAR        =fileNetCDF(istart:iend)


      s = nf90_create(fileNetCDF, or(nf90_clobber,NF90_HDF5), nc)

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
       !s = nf90_def_var(nc,'time'   , nf90_double,  (/timid/)  , idTim)

        s = nf90_def_var(nc,'TRB'//trim(VAR), nf90_double, (/xid,yid,depid,timid/), idB)
        s = nf90_def_var(nc,'TRN'//trim(VAR), nf90_double, (/xid,yid,depid,timid/), idN)

        s= nf90_put_att(nc,idTim ,'Units', 'seconds since 1582-10-15 00:00:00')
      
        s = nf90_put_att(nc,idB   , 'missing_value',1.e+20)
        s = nf90_put_att(nc,idN   , 'missing_value',1.e+20)
        s =nf90_enddef(nc)

        s = nf90_put_var(nc, idLon,  totglamt)
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idLat,  totgphit)
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idLev,     gdept)
       call handle_err1(s,counter,fileNetCDF)

        s = nf90_put_var(nc, idTim,    julian)
       call handle_err1(s,counter,fileNetCDF)

        s = nf90_put_var(nc, idB,      tottrb)
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idN,      tottrn)
       call handle_err1(s,counter,fileNetCDF)

        s =nf90_close(nc)


       END SUBROUTINE write_restart

       !****************************************************************************
       !****************************************************************************
       !****************************************************************************

       SUBROUTINE WRITE_AVE(fileNetCDF,VAR, datefrom, dateTo,M)
       USE netcdf
       USE myalloc
       IMPLICIT NONE

       CHARACTER*(*) fileNetCDF
       character(LEN=17) datefrom, dateTo
       real(4) M(jpk, jpjglo, jpiglo)

       character(LEN=20) VAR
       integer istart,iend

       integer s, nc, counter
       integer timid, depid, yid, xid
       integer idvartime,idgdept,idphit,idlamt,idVAR
       real lat_actual_range(2), lon_actual_range(2), depth_actual_range(2)
         lon_actual_range=(/-9.25  , 36.0   /)
         lat_actual_range=(/30.5   , 44.5   /)
       depth_actual_range=(/ 4.9991,4450.068/)

!      istart=index(fileNetCDF,'00.')+3
!      istart=index(fileNetCDF,'.nc')-3
!      iend  =index(fileNetCDF,'.nc')-1
!      VAR        =fileNetCDF(istart:iend)

        counter=0

        !s = nf90_create(fileNetCDF, or(nf90_clobber,NF90_HDF5), nc)
      
        s = nf90_create(fileNetCDF,or(or(nf90_clobber,NF90_HDF5),NF90_HDF5),nc)
        call handle_err1(s,counter,fileNetCDF)
        ! *********** GLOBAL ********************
        s = nf90_put_att(nc, nf90_global, 'Convenctions' ,'COARDS')
      call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_att(nc, nf90_global, 'DateStart'     , datefrom)
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_att(nc, nf90_global, 'Date__End'     ,   dateTo)
       call handle_err1(s,counter,fileNetCDF)

        ! *********** DIMENSIONS ****************
        s= nf90_def_dim(nc,'lon'           , jpiglo,  xid)
       call handle_err1(s,counter,fileNetCDF)
        s= nf90_def_dim(nc,'lat'           , jpjglo,  yid)
       call handle_err1(s,counter,fileNetCDF)
        s= nf90_def_dim(nc,'depth'         , jpk   ,depid)
       call handle_err1(s,counter,fileNetCDF)
        s= nf90_def_dim(nc,'time'  , NF90_UNLIMITED,timid)
       call handle_err1(s,counter,fileNetCDF)

        ! ********** VARIABLES *****************
        !!s = nf90_def_var(nc,'time',         nf90_double,(/timid/),       idvartime)
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_def_var(nc,'depth',        nf90_float, (/depid/),         idgdept)
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_def_var(nc,'lat'   ,       nf90_float, (/yid/),            idphit)
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_def_var(nc,'lon'   ,       nf90_float, (/xid/),            idlamt)
       call handle_err1(s,counter,fileNetCDF)
   
        s = nf90_def_var(nc,trim(VAR) ,        nf90_float, (/xid,yid,depid,timid/),  idVAR)
       call handle_err1(s,counter,fileNetCDF)

        s = nf90_put_att(nc,idgdept,'units'        ,'meter')
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_att(nc,idgdept,'positive'     ,'down')
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_att(nc,idgdept,'actual_range' ,depth_actual_range)
       call handle_err1(s,counter,fileNetCDF)

        s = nf90_put_att(nc,idphit, 'units'        ,'degrees_north')
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_att(nc,idphit, 'long_name'    ,'Latitude')
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_att(nc,idphit,'actual_range' ,lat_actual_range)
       call handle_err1(s,counter,fileNetCDF)

        s = nf90_put_att(nc,idlamt, 'units'        ,'degrees_east')
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_att(nc,idlamt, 'long_name'    ,'Longitude')
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_att(nc,idlamt,'actual_range' ,lon_actual_range)
       call handle_err1(s,counter,fileNetCDF)

        s = nf90_put_att(nc,idVAR, 'long_name'    ,VAR)
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_att(nc,idVAR, 'missing_value'    ,1.e+20)
       call handle_err1(s,counter,fileNetCDF)

        s =nf90_enddef(nc)
       call handle_err1(s,counter,fileNetCDF)


        s = nf90_put_var(nc, idlamt,   REAL(totglamt(jpjglo,:),4) )
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idphit,   REAL(totgphit(:,jpiglo),4) )
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idgdept,  REAL(   gdept,     4) )
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idVAR  ,  M )                    
       call handle_err1(s,counter,fileNetCDF)


        s =nf90_close(nc)
       call handle_err1(s,counter,fileNetCDF)

       END SUBROUTINE WRITE_AVE

       SUBROUTINE WRITE_AVE_2D(fileNetCDF,VAR, datefrom, dateTo,M)
       USE netcdf
       USE myalloc
       IMPLICIT NONE

       CHARACTER*(*) fileNetCDF
       character(LEN=17) datefrom, dateTo
       real(4) M(jpjglo, jpiglo)

       character(LEN=20) VAR
       integer istart,iend

       integer s, nc, counter
       integer timid, yid, xid
       integer idvartime,idphit,idlamt,idVAR
       real lat_actual_range(2), lon_actual_range(2)
         lon_actual_range=(/-9.25  , 36.0   /)
         lat_actual_range=(/30.5   , 44.5   /)


!      istart=index(fileNetCDF,'00.')+3
!      iend  =index(fileNetCDF,'.nc')-1
!      VAR        =fileNetCDF(istart:iend)


        s = nf90_create(fileNetCDF, or(nf90_clobber,NF90_HDF5), nc)
        ! *********** GLOBAL ********************
        s = nf90_put_att(nc, nf90_global, 'Convenctions'   ,'COARDS')
        s = nf90_put_att(nc, nf90_global, 'DateStart'     , datefrom)
        s = nf90_put_att(nc, nf90_global, 'Date__End'     ,   dateTo)

        ! *********** DIMENSIONS ****************
        s= nf90_def_dim(nc,'lon'           , jpiglo,  xid)
        s= nf90_def_dim(nc,'lat'           , jpjglo,  yid)
        s= nf90_def_dim(nc,'time'  , NF90_UNLIMITED,timid)

        ! ********** VARIABLES *****************
        !s = nf90_def_var(nc,'time',         nf90_double,(/timid/),       idvartime)
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
        ! epascolo warning
        s = nf90_put_var(nc, idlamt,   REAL(totglamt(jpjglo,:),4) )
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idphit,   REAL(totgphit(:,jpiglo),4) )
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idVAR  ,  M )                    
       call handle_err1(s,counter,fileNetCDF)


        s =nf90_close(nc)

       END SUBROUTINE WRITE_AVE_2D

!****************************************************************************
       !****************************************************************************
       !****************************************************************************

       SUBROUTINE WRITE_AVE_BKP(fileNetCDF, VAR,datefrom, dateTo,M, ave_counter)
       USE netcdf
       USE myalloc
       IMPLICIT NONE

       CHARACTER*(*) fileNetCDF
       character(LEN=17) datefrom, dateTo
       real(8) M(jpk, jpjglo, jpiglo)
       integer ave_counter

       !local
       character(LEN=20) VAR
       integer istart,iend

       integer s, nc, counter
       integer timid, depid, yid, xid
       integer idvartime,idgdept,idphit,idlamt,idVAR


!      istart=index(fileNetCDF,'00.')+3
!      iend  =index(fileNetCDF,'.nc')-1
!      VAR        =fileNetCDF(istart:iend)


        s = nf90_create(fileNetCDF, or(nf90_clobber,NF90_HDF5), nc)
        ! *********** GLOBAL ********************
        s = nf90_put_att(nc, nf90_global, 'Convenctions'  ,    'COARDS')
        s = nf90_put_att(nc, nf90_global, 'DateStart'     ,    datefrom)
        s = nf90_put_att(nc, nf90_global, 'Date__End'     ,      dateTo)
        s = nf90_put_att(nc, nf90_global, 'ave_counter'   , ave_counter)

        ! *********** DIMENSIONS ****************
        s= nf90_def_dim(nc,'lon'           , jpiglo,  xid)
        s= nf90_def_dim(nc,'lat'           , jpjglo,  yid)
        s= nf90_def_dim(nc,'depth'         , jpk   ,depid)
        s= nf90_def_dim(nc,'time'  , NF90_UNLIMITED,timid)

        ! ********** VARIABLES *****************
        !s = nf90_def_var(nc,'time',         nf90_double,(/timid/),       idvartime)
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

        ! epascolo warning
        s = nf90_put_var(nc, idlamt,   REAL(totglamt(jpjglo,:),4) )
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idphit,   REAL(totgphit(:,jpiglo),4) )
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idgdept,  REAL(   gdept,          4) )
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idVAR  ,                          M  )
       call handle_err1(s,counter,fileNetCDF)


        s =nf90_close(nc)

       END SUBROUTINE WRITE_AVE_BKP

      SUBROUTINE WRITE_AVE_2d_BKP(fileNetCDF,VAR, datefrom, dateTo,M, ave_counter)
       USE netcdf
       USE myalloc
       IMPLICIT NONE

       CHARACTER*(*) fileNetCDF
       character(LEN=17) datefrom, dateTo
       real(8) M(jpjglo, jpiglo)
       integer ave_counter

       !local
       character(LEN=20) VAR
       integer istart,iend

       integer s, nc, counter
       integer timid, yid, xid
       integer idvartime,idgdept,idphit,idlamt,idVAR


!      istart=index(fileNetCDF,'00.')+3
!      iend  =index(fileNetCDF,'.nc')-1
!      VAR        =fileNetCDF(istart:iend)


        s = nf90_create(fileNetCDF, or(nf90_clobber,NF90_HDF5), nc)
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
        !s = nf90_def_var(nc,'time',         nf90_double,(/timid/),       idvartime)
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


        s = nf90_put_var(nc, idlamt,   REAL(totglamt(jpjglo,:),4) )
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idphit,   REAL(totgphit(:,jpiglo),4) )
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idVAR  ,                          M  )
       call handle_err1(s,counter,fileNetCDF)


        s =nf90_close(nc)

       END SUBROUTINE WRITE_AVE_2d_BKP



       !****************************************************************************
       ! physical data were printed as grid3D.dat grid2D.dat files
       SUBROUTINE physDump(fileNetCDF, datefrom, dateTo)
       USE myalloc
       USE netcdf
       USE IO_mem, only: d2f3d, d2f2d

       IMPLICIT NONE

        character fileNetCDF*(*)
        character(LEN=17) datefrom, dateTo
        integer s, nc, counter
        integer timid, depid, yid, xid
        integer idvartime, idgdept, idphit, idlamt
        integer idT, idS, idU, idV, idW, idEddy,ide3t, idR, idWs, idE


        s = nf90_create(fileNetCDF, or(nf90_clobber,NF90_HDF5), nc)
        ! *********** GLOBAL ********************
        s = nf90_put_att(nc, nf90_global, 'Convenctions'   ,'COARDS')
        s = nf90_put_att(nc, nf90_global, 'DateStart'     , datefrom)
        s = nf90_put_att(nc, nf90_global, 'Date__End'     ,   dateTo)

        ! *********** DIMENSIONS ****************
        s= nf90_def_dim(nc,'x'           , jpiglo,  xid)
        s= nf90_def_dim(nc,'y'           , jpjglo,  yid)
        s= nf90_def_dim(nc,'deptht'      , jpk   ,depid)
        s= nf90_def_dim(nc,'time_counter', NF90_UNLIMITED,timid)


        ! ********** VARIABLES *****************
        !s = nf90_def_var(nc,'time_counter',  nf90_double,(/timid/),           idvartime)
        s = nf90_def_var(nc,'deptht',        nf90_float, (/depid/),             idgdept)
        s = nf90_def_var(nc,'nav_lat',       nf90_float, (/yid/),                idphit)
        s = nf90_def_var(nc,'nav_lon',       nf90_float, (/xid/),                idlamt)

        s = nf90_def_var(nc,'vosaline' ,     nf90_float, (/xid,yid,depid,timid/),   idS)
        s = nf90_def_var(nc,'votemper' ,     nf90_float, (/xid,yid,depid,timid/),   idT)
        s = nf90_def_var(nc,'vozocrtx' ,     nf90_float, (/xid,yid,depid,timid/),   idU)
        s = nf90_def_var(nc,'vomecrty' ,     nf90_float, (/xid,yid,depid,timid/),   idV)
        s = nf90_def_var(nc,'vovecrtz' ,     nf90_float, (/xid,yid,depid,timid/),   idW)
        s = nf90_def_var(nc,'votkeavt' ,     nf90_float, (/xid,yid,depid,timid/),idEddy)
        s = nf90_def_var(nc,'e3t'      ,     nf90_double, (/xid,yid,depid,timid/), ide3t)

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

        s = nf90_put_att(nc,ide3t , 'long_name'    ,'Vartical scale factor')
        s = nf90_put_att(nc,ide3t , 'units'        ,'m')
        s = nf90_put_att(nc,ide3t , 'missing_value',1.e+20)

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

        s = nf90_put_var(nc, idlamt,  REAL(totglamt(jpjglo,:),4)  )
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idphit,  REAL(totgphit(:,jpiglo),4)  )
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idgdept, REAL(   gdept          ,4)  )
       call handle_err1(s,counter,fileNetCDF)
!      3D vars
        d2f3d = REAL(totsnIO,4) 
       s = nf90_put_var(nc, idS,    d2f3d) 
        call handle_err1(s,counter,fileNetCDF)
        d2f3d = REAL(tottnIO,4) 
       s = nf90_put_var(nc, idT,    d2f3d) 
        call handle_err1(s,counter,fileNetCDF)
        d2f3d = REAL(totunIO,4) 
       s = nf90_put_var(nc, idU,    d2f3d) 
        call handle_err1(s,counter,fileNetCDF)
        d2f3d = REAL(totvnIO,4) 
       s = nf90_put_var(nc, idV,    d2f3d) 
        call handle_err1(s,counter,fileNetCDF)
        d2f3d = REAL(totwnIO,4) 
       s = nf90_put_var(nc, idW,    d2f3d) 
        call handle_err1(s,counter,fileNetCDF)
        d2f3d = REAL(totavtIO,4)
       s = nf90_put_var(nc, idEddy, d2f3d) 
        call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, ide3t, tote3tIO) 
        call handle_err1(s,counter,fileNetCDF)
!       d2f3d = REAL(tote3tIO,4)
       s = nf90_put_var(nc, ide3t, d2f3d) 
        call handle_err1(s,counter,fileNetCDF)

!       2D vars
        d2f2d =REAL(totvatmIO,4)
       s = nf90_put_var(nc, idWs,   d2f2d) 
       call handle_err1(s,counter,fileNetCDF)
        d2f2d =REAL(totempIO ,4)
       s = nf90_put_var(nc, idE,    d2f2d) 
       call handle_err1(s,counter,fileNetCDF)
        d2f2d =REAL(totqsrIO ,4)
       s = nf90_put_var(nc, idR,    d2f2d) 
       call handle_err1(s,counter,fileNetCDF)



        s= nf90_close(nc)


       END SUBROUTINE physDump


       !****************************************************************************
       ! physical data were printed as grid3D.dat grid2D.dat files
       SUBROUTINE physDump_bkp(fileNetCDF, datefrom, dateTo,ave_counter)
       USE myalloc
       USE netcdf


       IMPLICIT NONE

        character fileNetCDF*(*)
        character(LEN=17) datefrom, dateTo
        integer ave_counter
        integer s, nc, counter
        integer timid, depid, yid, xid
        integer idvartime, idgdept, idphit, idlamt
        integer idT, idS, idU, idV, idW, idEddy, ide3t, idR, idWs, idE


        s = nf90_create(fileNetCDF, or(nf90_clobber,NF90_HDF5), nc)
        ! *********** GLOBAL ********************
        s = nf90_put_att(nc, nf90_global, 'Convenctions'   ,   'COARDS')
        s = nf90_put_att(nc, nf90_global, 'DateStart'     ,    datefrom)
        s = nf90_put_att(nc, nf90_global, 'Date__End'     ,      dateTo)
        s = nf90_put_att(nc, nf90_global, 'ave_counter'   , ave_counter)

        ! *********** DIMENSIONS ****************
        s= nf90_def_dim(nc,'x'           , jpiglo,  xid)
        s= nf90_def_dim(nc,'y'           , jpjglo,  yid)
        s= nf90_def_dim(nc,'deptht'      , jpk   ,depid)
        s= nf90_def_dim(nc,'time_counter', NF90_UNLIMITED,timid)


        ! ********** VARIABLES *****************
        s = nf90_def_var(nc,'time_counter',  nf90_double,(/timid/),           idvartime)
        s = nf90_def_var(nc,'deptht',        nf90_float, (/depid/),             idgdept)
        s = nf90_def_var(nc,'nav_lat',       nf90_float, (/yid/),                idphit)
        s = nf90_def_var(nc,'nav_lon',       nf90_float, (/xid/),                idlamt)

        s = nf90_def_var(nc,'vosaline' ,     nf90_double, (/xid,yid,depid,timid/),   idS)
        s = nf90_def_var(nc,'votemper' ,     nf90_double, (/xid,yid,depid,timid/),   idT)
        s = nf90_def_var(nc,'vozocrtx' ,     nf90_double, (/xid,yid,depid,timid/),   idU)
        s = nf90_def_var(nc,'vomecrty' ,     nf90_double, (/xid,yid,depid,timid/),   idV)
        s = nf90_def_var(nc,'vovecrtz' ,     nf90_double, (/xid,yid,depid,timid/),   idW)
        s = nf90_def_var(nc,'votkeavt' ,     nf90_double, (/xid,yid,depid,timid/),idEddy)
        s = nf90_def_var(nc,'e3t'      ,     nf90_double, (/xid,yid,depid,timid/), ide3t)

        s = nf90_def_var(nc,'soshfldo' ,     nf90_double, (/xid,yid,      timid/),   idR)
        s = nf90_def_var(nc,'sowindsp' ,     nf90_double, (/xid,yid,      timid/),  idWs)
        s = nf90_def_var(nc,'sowaflcd' ,     nf90_double, (/xid,yid,      timid/),   idE)


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

        s = nf90_put_att(nc,ide3t , 'long_name'    ,'Vertical scale factor')
        s = nf90_put_att(nc,ide3t , 'units'        ,'m')
        s = nf90_put_att(nc,ide3t , 'missing_value',1.e+20)

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
       call handle_err1(s,counter,fileNetCDF)

        counter=0

        s = nf90_put_var(nc, idlamt,  REAL(totglamt(jpjglo,:),4)  )
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idphit,  REAL(totgphit(:,jpiglo),4)  )
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idgdept, REAL(   gdept          ,4)  )
       call handle_err1(s,counter,fileNetCDF)
!      3D vars
        s = nf90_put_var(nc, idS,    totsnIO) 
        call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idT,    tottnIO) 
        call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idU,    totunIO) 
        call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idV,    totvnIO) 
        call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idW,    totwnIO) 
        call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idEddy,totavtIO) 
        call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, ide3t ,tote3tIO) 
        call handle_err1(s,counter,fileNetCDF)

!       2D vars
        s = nf90_put_var(nc, idWs,  totvatmIO) 
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idE,    totempIO) 
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idR,    totqsrIO) 
       call handle_err1(s,counter,fileNetCDF)



        s= nf90_close(nc)


       END SUBROUTINE physDump_bkp





       !****************************************************************************
       !****************************************************************************


        SUBROUTINE get_att_char(fileNetCDF,attname,attvalue)
        use netcdf
        IMPLICIT NONE
        character fileNetCDF*(*)
        character attname*(*)
        character attvalue*(*)

        integer stat, ncid
        integer counter
        counter=0
        stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)  
             call handle_err1(stat, counter,FileNetCDF)
        stat= nf90_get_att(ncid, nf90_global, attname, attvalue)
       call handle_err1(stat, counter,FileNetCDF)
        stat = nf90_close(ncid)
      

        END SUBROUTINE get_att_char



        !****************************************************************************
        !****************************************************************************

        SUBROUTINE get_att_int(fileNetCDF,attname,attvalue)
        use netcdf
        IMPLICIT NONE
        character fileNetCDF*(*)
        character attname*(*)
        integer attvalue

        integer stat, ncid
        integer counter
        counter=0
        stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)  
             call handle_err1(stat, counter,FileNetCDF)
        stat= nf90_get_att(ncid, nf90_global, attname, attvalue)
       call handle_err1(stat, counter,FileNetCDF)
        stat = nf90_close(ncid)
      

        END SUBROUTINE get_att_int


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
