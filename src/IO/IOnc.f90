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
    SUBROUTINE readnc_scalar_double(fileNetCDF,varname, M)
    USE netcdf
    implicit none

    character,intent(in) :: fileNetCDF*(*) ,varname*(*)
    double precision,intent(inout) ::  M
    integer ncid, stat, VARid
    integer counter

    counter=0
    stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)
    call handle_err1(stat, counter,FileNetCDF)

    stat = nf90_inq_varid (ncid, varname, VARid)
    call handle_err2(stat, fileNetCDF,varname)
    call handle_err1(stat, counter,FileNetCDF)
    stat = nf90_get_var (ncid,VARid,M)
    stat = nf90_close(ncid)
    call handle_err1(stat, counter,FileNetCDF)

    END SUBROUTINE readnc_scalar_double
! ! *****************************************************************************

      SUBROUTINE readnc_slice_double(fileNetCDF,varname, M)
      USE myalloc
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
       call handle_err2(stat, fileNetCDF,varname)
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
      SUBROUTINE readnc_slice_int1(fileNetCDF,varname, M)
      USE myalloc
      USE netcdf
      implicit none


      character,intent(in) :: fileNetCDF*(*) ,varname*(*)
      INTEGER(kind=1),intent(inout) ::  M(jpk,jpj,jpi)
      
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
       call handle_err2(stat, fileNetCDF,varname)
       call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_get_var (ncid,VARid,copy_in,start, thecount)
      
      call handle_err2(stat, fileNetCDF,varname)        
      call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_close(ncid)                           
      call handle_err1(stat, counter,FileNetCDF)

      DO i=1,jpi
        DO j=1,jpj
          DO k=1,jpk
            M(k,j,i) = INT(copy_in(i,j,k),1)
          ENDDO
        ENDDO
      ENDDO   
      
      deallocate(copy_in)

      END SUBROUTINE readnc_slice_int1

! ********************************************************************************************

      subroutine readnc_slice_logical(fileNetCDF, varname, M)

          use myalloc
          use netcdf

          implicit none

          character, intent(in) :: fileNetCDF*(*), varname*(*)
          integer(1), intent(inout) :: M(jpk, jpj, jpi)

          integer(1), allocatable, dimension(:, :, :) :: copy_in
          integer ncid, stat, VARid, i, j, k
          integer counter
          integer thecount(4), start(4)

          allocate(copy_in(jpi, jpj, jpk))
          counter = 0
          start = (/ nimpp, njmpp, 1, 1 /)
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

          do i = 1, jpi
              do j = 1, jpj
                  do k = 1, jpk
                      M(k, j, i) = copy_in(i, j, k)
                  enddo
              enddo
          enddo

          deallocate(copy_in)

      end subroutine readnc_slice_logical

! ********************************************************************************************

      SUBROUTINE readnc_slice_float(fileNetCDF,varname, M, shift)
      USE myalloc
      USE netcdf
      implicit none


      character,intent(in) :: fileNetCDF*(*) ,varname*(*)
      integer, intent(in)  :: shift
      double precision,intent(inout) ::  M(jpk,jpj,jpi)
      
      real,allocatable,dimension(:,:,:) :: copy_in
      integer ncid, stat, VARid,i,j,k
      integer counter
      integer thecount(4), start(4)

      allocate(copy_in(jpi,jpj,jpk))
      counter = 0
      start    = (/nimpp+shift, njmpp,  1,  1/)
      thecount = (/jpi,           jpj, jpk, 1/)

      stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)  
       call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_inq_varid (ncid, varname, VARid)
       call handle_err2(stat, fileNetCDF,varname)
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
      integer,intent(inout),dimension(jpk,jpj,jpi) ::  M
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
      call handle_err2(stat, fileNetCDF,varname)
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
       call handle_err2(stat, fileNetCDF,varname)
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

     SUBROUTINE readnc_slice_float_2d(fileNetCDF,varname, M,shift)
     USE myalloc
      USE netcdf
      implicit none

      character ,intent(in) :: fileNetCDF*(*) ,varname*(*)
      integer , intent(in)  :: shift
      double precision,intent(inout) :: M(jpj,jpi)
      real,allocatable,dimension(:,:) :: copy_in
      
      integer ncid, stat, VARid,i,j,k
      integer counter
      integer thecount(3), start(3)


      allocate(copy_in(jpi,jpj))
      counter = 0
      start    = (/nimpp+shift, njmpp,  1/)
      thecount = (/jpi      ,     jpj,  1/)


      stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)  
       call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_inq_varid (ncid, varname, VARid)
       call handle_err2(stat, fileNetCDF,varname)
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
       call handle_err2(stat, fileNetCDF,varname)
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
      thecount = (/jpiglo,  jpjglo, 1, 1/)

      stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)  
      call handle_err1(stat, counter,FileNetCDF)
      stat = nf90_inq_varid (ncid, varname, VARid)
      call handle_err2(stat, fileNetCDF,varname)
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
       call handle_err2(stat, fileNetCDF,varname)
       call handle_err1(stat, counter,FileNetCDF)
        stat = nf90_get_var (ncid,VARid,ARRAY,start, thecount)
      
        call handle_err2(stat, fileNetCDF,varname)       
       call handle_err1(stat, counter,FileNetCDF)
        stat = nf90_close(ncid)                          
       call handle_err1(stat, counter,FileNetCDF)


        end SUBROUTINE readmask_double_1d


        !****************************************************************************
        SUBROUTINE readnc_double_1d(fileNetCDF,varname,dim1,ARRAY)
        use netcdf
        use myalloc
        implicit none

        character,intent(in) :: fileNetCDF*(*) ,varname*(*)
        integer,intent(in) :: dim1
        double precision,intent(inout),dimension(dim1) :: ARRAY
        
        integer ncid, stat, VARid
        integer counter
        
        counter=0

        stat = nf90_open(fileNetCDF, nf90_nowrite, ncid) 
       call handle_err1(stat, counter,FileNetCDF)
        stat = nf90_inq_varid (ncid, varname, VARid)
       call handle_err2(stat, fileNetCDF,varname)
       call handle_err1(stat, counter,FileNetCDF)
        stat = nf90_get_var (ncid,VARid,ARRAY)           
      
        call handle_err2(stat, fileNetCDF,varname)       
       call handle_err1(stat, counter,FileNetCDF)
        stat = nf90_close(ncid)                          
       call handle_err1(stat, counter,FileNetCDF)


        end SUBROUTINE readnc_double_1d

        !****************************************************************************
        SUBROUTINE readnc_int_1d(fileNetCDF,varname,dim1,ARRAY)
        use netcdf
        use myalloc
        implicit none

        character,intent(in) :: fileNetCDF*(*) ,varname*(*)
        integer,intent(in) :: dim1
        integer,intent(inout),dimension(dim1) :: ARRAY

        integer ncid, stat, VARid
        integer counter
        

        counter=0

       stat = nf90_open(fileNetCDF, nf90_nowrite, ncid) 
       call handle_err1(stat, counter,FileNetCDF)
       stat = nf90_inq_varid (ncid, varname, VARid)
       call handle_err2(stat, fileNetCDF,varname)
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

       SUBROUTINE write_restart(fileNetCDF,VAR, julian,deflate, deflate_level)
       USE netcdf
       USE myalloc

       IMPLICIT NONE
       CHARACTER*(*),intent(in) :: fileNetCDF
       double precision,intent(in) :: julian
       CHARACTER(*),intent(in) ::  VAR
       integer, intent(in) :: deflate, deflate_level

       ! local
       CHARACTER(LEN=17) :: TimeString
       integer :: istart, iend
       integer :: s, nc, counter
       integer :: timid, depid, yid, xid, xaid, yaid, zaid
       integer :: idB, idN, idLon, idLat, idLev, idTim
       integer shuffle
       double precision,allocatable,dimension(:,:,:) :: copy_in
       TimeString =fileNetCDF(14:30)
       shuffle       = 0

       ! Just to try without 'or'
       ! s = nf90_create(fileNetCDF, or(nf90_clobber,NF90_HDF5), nc)
       s = nf90_create(fileNetCDF, NF90_HDF5, nc)

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
        s = nf90_def_var(nc,'TRN'//VAR, nf90_double, (/xid,yid,depid,timid/), idN)
        s = nf90_def_var_deflate(nc, idN, shuffle, deflate, deflate_level)
        call handle_err1(s,counter,fileNetCDF)
        !s= nf90_put_att(nc,idTim ,'Units', 'seconds since 1582-10-15 00:00:00')
      
        s = nf90_put_att(nc,idN   , 'missing_value',1.e+20)
        s =nf90_enddef(nc)
        s = nf90_put_var(nc, idLon,  TRANSPOSE(totglamt))
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idLat,  TRANSPOSE(totgphit))
       call handle_err1(s,counter,fileNetCDF)
       

        s = nf90_put_var(nc, idLev,     gdept)
       call handle_err1(s,counter,fileNetCDF)

       allocate(copy_in(jpiglo, jpjglo, jpk))

       call switch_index_double(tottrn,copy_in,jpiglo,jpjglo,jpk)
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idN,      copy_in)
       call handle_err1(s,counter,fileNetCDF)
       deallocate(copy_in)
        s =nf90_close(nc)


       END SUBROUTINE write_restart

       !****************************************************************************
       !****************************************************************************
       !****************************************************************************

       SUBROUTINE WRITE_AVE(fileNetCDF,VAR, datefrom, dateTo,M,deflate, deflate_level)
       USE netcdf
       USE myalloc
       IMPLICIT NONE

       CHARACTER*(*),intent(in) :: fileNetCDF
       character(LEN=20),intent(in) :: VAR
       character(LEN=17),intent(in) :: datefrom, dateTo
       double precision, dimension(jpk, jpjglo, jpiglo),intent(in) :: M
       integer, intent(in) :: deflate, deflate_level
       
       real,allocatable,dimension(:,:,:) :: copy_in
       integer istart,iend
       integer s, nc, counter
       integer timid, depid, yid, xid
       integer idvartime,idgdept,idphit,idlamt,idVAR
       integer shuffle
       real lat_actual_range(2), lon_actual_range(2), depth_actual_range(2)
         lon_actual_range=(/-9.25  , 36.0   /)
         lat_actual_range=(/30.5   , 44.5   /)
       depth_actual_range=(/ 4.9991,4450.068/)
       shuffle       = 0


        counter=0

        !s = nf90_create(fileNetCDF, or(nf90_clobber,NF90_HDF5), nc)
      
        ! Just to try withour 'or'
        ! s = nf90_create(fileNetCDF,or(or(nf90_clobber,NF90_HDF5),NF90_HDF5),nc)
        s = nf90_create(fileNetCDF, NF90_HDF5, nc)
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
        !call handle_err1(s,counter,fileNetCDF)
        s = nf90_def_var(nc,'depth',        nf90_float, (/depid/),         idgdept)
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_def_var(nc,'lat'   ,       nf90_float, (/yid/),            idphit)
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_def_var(nc,'lon'   ,       nf90_float, (/xid/),            idlamt)
       call handle_err1(s,counter,fileNetCDF)
   
        s = nf90_def_var(nc,trim(VAR) ,        nf90_float, (/xid,yid,depid,timid/),  idVAR)
       call handle_err1(s,counter,fileNetCDF)
       s = nf90_def_var_deflate(nc, idVAR, shuffle, deflate, deflate_level)
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

       allocate(copy_in(jpiglo, jpjglo, jpk))
       call switch_index_rout(M,copy_in,jpiglo,jpjglo,jpk)
        s = nf90_put_var(nc, idVAR  ,  copy_in )                    
       call handle_err1(s,counter,fileNetCDF)
       deallocate(copy_in)

        s =nf90_close(nc)
       call handle_err1(s,counter,fileNetCDF)

       END SUBROUTINE WRITE_AVE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       SUBROUTINE WRITE_AVE_2D(fileNetCDF,VAR, datefrom, dateTo,M)
       USE netcdf
       USE myalloc
       IMPLICIT NONE

       character*(*),intent(in) :: fileNetCDF
       character(LEN=17),intent(in) :: datefrom, dateTo
       real,intent(in),dimension(jpjglo, jpiglo) :: M

       character(LEN=20) :: VAR
       integer :: istart,iend

       integer :: s, nc, counter
       integer :: timid, yid, xid
       integer :: idvartime,idphit,idlamt,idVAR
       real :: lat_actual_range(2), lon_actual_range(2)
         lon_actual_range=(/-9.25  , 36.0   /)
         lat_actual_range=(/30.5   , 44.5   /)



        ! Just to try without 'or'
        ! s = nf90_create(fileNetCDF, or(nf90_clobber,NF90_HDF5), nc)
        s = nf90_create(fileNetCDF, NF90_HDF5, nc)
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
        s = nf90_put_var(nc, idlamt,  REAL(totglamt(jpjglo,:),4) )
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idphit,  REAL(totgphit(:,jpiglo),4) )
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idVAR  ,  transpose(M) )                    
       call handle_err1(s,counter,fileNetCDF)


        s =nf90_close(nc)

       END SUBROUTINE WRITE_AVE_2D

!****************************************************************************
       !****************************************************************************
       !****************************************************************************

       SUBROUTINE WRITE_AVE_BKP(fileNetCDF, VAR,datefrom, dateTo,M, elapsed_time, deflate, deflate_level)
       USE netcdf
       USE myalloc
       IMPLICIT NONE

       CHARACTER*(*),intent(in) :: fileNetCDF
       character(LEN=20), intent(in):: VAR
       character(LEN=17),intent(in) :: datefrom, dateTo
       double precision,dimension(jpk, jpjglo, jpiglo),intent(in) :: M
       double precision,intent(in) :: elapsed_time
       integer, intent(in) :: deflate, deflate_level


       !local
       double precision, allocatable,dimension(:,:,:) :: copy_in
       integer istart,iend

       integer s, nc, counter
       integer timid, depid, yid, xid
       integer idvartime,idgdept,idphit,idlamt,idVAR
       integer shuffle
       shuffle       = 0



        ! Just to try without 'or'
        ! s = nf90_create(fileNetCDF, or(nf90_clobber,NF90_HDF5), nc)
        s = nf90_create(fileNetCDF, NF90_HDF5, nc)
        ! *********** GLOBAL ********************
        s = nf90_put_att(nc, nf90_global, 'Convenctions'  ,    'COARDS')
        s = nf90_put_att(nc, nf90_global, 'DateStart'     ,    datefrom)
        s = nf90_put_att(nc, nf90_global, 'Date__End'     ,      dateTo)


        ! *********** DIMENSIONS ****************
        s= nf90_def_dim(nc,'lon'           , jpiglo,  xid)
        s= nf90_def_dim(nc,'lat'           , jpjglo,  yid)
        s= nf90_def_dim(nc,'depth'         , jpk   ,depid)
        s= nf90_def_dim(nc,'time'          , 1     ,timid)

        ! ********** VARIABLES *****************
        s = nf90_def_var(nc,'elapsed_time', nf90_double,(/timid/),       idvartime)
        s = nf90_def_var(nc,'depth',        nf90_float, (/depid/),         idgdept)
        s = nf90_def_var(nc,'lat'   ,       nf90_float, (/yid/),            idphit)
        s = nf90_def_var(nc,'lon'   ,       nf90_float, (/xid/),            idlamt)

       s = nf90_def_var(nc,trim(VAR) , nf90_double,(/xid,yid,depid,timid/),  idVAR)
       s = nf90_def_var_deflate(nc, idVAR, shuffle, deflate, deflate_level)
       call handle_err1(s,counter,fileNetCDF)

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

        s = nf90_put_var(nc,idvartime, elapsed_time)
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idlamt,   REAL(totglamt(jpjglo,:),4) )
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idphit,   REAL(totgphit(:,jpiglo),4) )
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idgdept,  REAL(   gdept,          4) )
       call handle_err1(s,counter,fileNetCDF)

       allocate(copy_in(jpiglo, jpjglo, jpk))
       call switch_index_double(M,copy_in,jpiglo,jpjglo,jpk)
       s = nf90_put_var(nc, idVAR  , copy_in  )
       deallocate(copy_in)

       call handle_err1(s,counter,fileNetCDF)


        s =nf90_close(nc)

       END SUBROUTINE WRITE_AVE_BKP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      SUBROUTINE WRITE_AVE_2d_BKP(fileNetCDF,VAR, datefrom, dateTo,M, elapsed_time)
       USE netcdf
       USE myalloc
       IMPLICIT NONE

       CHARACTER*(*),intent(in) :: fileNetCDF
       character(LEN=17),intent(in) :: datefrom, dateTo
       double precision,intent(in),dimension(jpjglo, jpiglo) :: M
       double precision,intent(in) :: elapsed_time

       !local
       character(LEN=20) VAR
       integer istart,iend

       integer s, nc, counter
       integer timid, yid, xid
       integer idvartime,idgdept,idphit,idlamt,idVAR 




        ! Just to try without 'or'
        ! s = nf90_create(fileNetCDF, or(nf90_clobber,NF90_HDF5), nc)
        s = nf90_create(fileNetCDF, NF90_HDF5, nc)
        ! *********** GLOBAL ********************
        s = nf90_put_att(nc, nf90_global, 'Convenctions'  ,    'COARDS')
        s = nf90_put_att(nc, nf90_global, 'DateStart'     ,    datefrom)
        s = nf90_put_att(nc, nf90_global, 'Date__End'     ,      dateTo)

        ! *********** DIMENSIONS ****************
        s= nf90_def_dim(nc,'lon'           , jpiglo,  xid)
        s= nf90_def_dim(nc,'lat'           , jpjglo,  yid)
        s= nf90_def_dim(nc,'time'          , 1     ,timid)

        ! ********** VARIABLES *****************
        s = nf90_def_var(nc,'elapsed_time', nf90_double,(/timid/),       idvartime)
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

        s = nf90_put_var(nc,idvartime, elapsed_time)
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idlamt,   REAL(totglamt(jpjglo,:),4) )
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idphit,   REAL(totgphit(:,jpiglo),4) )
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idVAR, transpose(M))
       call handle_err1(s,counter,fileNetCDF)


        s =nf90_close(nc)

       END SUBROUTINE WRITE_AVE_2d_BKP



       !****************************************************************************
       SUBROUTINE physDump(fileNetCDF, datefrom, dateTo)
       USE myalloc
       USE netcdf
       USE IO_mem, only: d2f2d

       IMPLICIT NONE

        character fileNetCDF*(*)
        character(LEN=17) datefrom, dateTo
        integer s, nc, counter
        integer timid, depid, yid, xid
        integer idvartime, idgdept, idphit, idlamt
        integer idT, idS, idU, idV, idW, idEddy,ide3t, idR, idWs, idE
        real,allocatable,dimension(:,:,:) :: copy_in
        real,allocatable,dimension(:,:) :: copy_in_2d
        allocate(copy_in(jpiglo,jpjglo,jpk))
        allocate(copy_in_2d(jpiglo,jpjglo))

        ! Just to try without 'or'
        ! s = nf90_create(fileNetCDF, or(nf90_clobber,NF90_HDF5), nc)
        s = nf90_create(fileNetCDF, NF90_HDF5, nc)
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
        
        call switch_index_rout(totsnIO,copy_in,jpiglo,jpjglo,jpk)
        s = nf90_put_var(nc, idS,    copy_in) 
        call handle_err1(s,counter,fileNetCDF)

      
        call switch_index_rout(tottnIO,copy_in,jpiglo,jpjglo,jpk)
       s = nf90_put_var(nc, idT,    copy_in) 
        call handle_err1(s,counter,fileNetCDF)

     
        call switch_index_rout(totunIO,copy_in,jpiglo,jpjglo,jpk) 
       s = nf90_put_var(nc, idU,    copy_in) 
        call handle_err1(s,counter,fileNetCDF)

        call switch_index_rout(totvnIO,copy_in,jpiglo,jpjglo,jpk) 
       s = nf90_put_var(nc, idV,    copy_in) 
        call handle_err1(s,counter,fileNetCDF)

     
        call switch_index_rout(totwnIO,copy_in,jpiglo,jpjglo,jpk) 
       s = nf90_put_var(nc, idW,    copy_in) 
        call handle_err1(s,counter,fileNetCDF)

        call switch_index_rout(totavtIO,copy_in,jpiglo,jpjglo,jpk)
       s = nf90_put_var(nc, idEddy, copy_in) 
        call handle_err1(s,counter,fileNetCDF)

        call switch_index_rout(tote3tIO,copy_in,jpiglo,jpjglo,jpk)
        s = nf90_put_var(nc, ide3t, copy_in) 
        call handle_err1(s,counter,fileNetCDF)
!       2D vars
        copy_in_2d =transpose(REAL(totvatmIO,4))
       s = nf90_put_var(nc, idWs,   copy_in_2d) 
       call handle_err1(s,counter,fileNetCDF)
        copy_in_2d =transpose(REAL(totempIO ,4))
       s = nf90_put_var(nc, idE,    copy_in_2d) 
       call handle_err1(s,counter,fileNetCDF)
        copy_in_2d =transpose(REAL(totqsrIO ,4))
       s = nf90_put_var(nc, idR,    copy_in_2d) 
       call handle_err1(s,counter,fileNetCDF)



        s= nf90_close(nc)
        deallocate(copy_in_2d)
        deallocate(copy_in)


       END SUBROUTINE physDump


       !****************************************************************************
       ! physical data were printed as grid3D.dat grid2D.dat files
       SUBROUTINE physDump_bkp(fileNetCDF, datefrom, dateTo,elapsed_time)
       USE myalloc
       USE netcdf


       IMPLICIT NONE

        character fileNetCDF*(*)
        character(LEN=17) datefrom, dateTo
        double precision,intent(in) :: elapsed_time
        integer s, nc, counter
        integer timid, depid, yid, xid
        integer  idvartime,idgdept, idphit, idlamt
        integer idT, idS, idU, idV, idW, idEddy, ide3t, idR, idWs, idE
        double precision,allocatable,dimension(:,:,:) :: copy_in
        double precision,allocatable,dimension(:,:) :: copy_in_2d
        allocate(copy_in(jpiglo,jpjglo,jpk))
        allocate(copy_in_2d(jpiglo,jpjglo))

        ! Just to try without 'or'
        ! s = nf90_create(fileNetCDF, or(nf90_clobber,NF90_HDF5), nc)
        s = nf90_create(fileNetCDF, NF90_HDF5, nc)
        ! *********** GLOBAL ********************
        s = nf90_put_att(nc, nf90_global, 'Convenctions'   ,   'COARDS')
        s = nf90_put_att(nc, nf90_global, 'DateStart'     ,    datefrom)
        s = nf90_put_att(nc, nf90_global, 'Date__End'     ,      dateTo)

        ! *********** DIMENSIONS ****************
        s= nf90_def_dim(nc,'x'           , jpiglo,  xid)
        s= nf90_def_dim(nc,'y'           , jpjglo,  yid)
        s= nf90_def_dim(nc,'deptht'      , jpk   ,depid)
        s= nf90_def_dim(nc,'time_counter', 1     ,timid)

        
        ! ********** VARIABLES *****************
        s = nf90_def_var(nc,'elapsed_time',  nf90_double,(/timid/),           idvartime)
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

        s = nf90_put_var(nc,idvartime, elapsed_time)
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idlamt,  REAL(totglamt(jpjglo,:),4)  )
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idphit,  REAL(totgphit(:,jpiglo),4)  )
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idgdept, REAL(   gdept          ,4)  )
       call handle_err1(s,counter,fileNetCDF)

!      3D vars
        call switch_index_double(totsnIO,copy_in,jpiglo,jpjglo,jpk)
        s = nf90_put_var(nc, idS,    copy_in) 
        call handle_err1(s,counter,fileNetCDF)

         call switch_index_double(tottnIO,copy_in,jpiglo,jpjglo,jpk)
        s = nf90_put_var(nc, idT,    copy_in) 
        call handle_err1(s,counter,fileNetCDF)

        call switch_index_double(totunIO,copy_in,jpiglo,jpjglo,jpk) 
        s = nf90_put_var(nc, idU,    copy_in) 
        call handle_err1(s,counter,fileNetCDF)

        call switch_index_double(totvnIO,copy_in,jpiglo,jpjglo,jpk)
        s = nf90_put_var(nc, idV,    copy_in) 
        call handle_err1(s,counter,fileNetCDF)

        call switch_index_double(totwnIO,copy_in,jpiglo,jpjglo,jpk)
        s = nf90_put_var(nc, idW,    copy_in) 
        call handle_err1(s,counter,fileNetCDF)

        call switch_index_double(totavtIO,copy_in,jpiglo,jpjglo,jpk)
        s = nf90_put_var(nc, idEddy,copy_in) 
        call handle_err1(s,counter,fileNetCDF)

        call switch_index_double(tote3tIO,copy_in,jpiglo,jpjglo,jpk)
        s = nf90_put_var(nc, ide3t ,copy_in) 
        call handle_err1(s,counter,fileNetCDF)

!       2D vars
        copy_in_2d =transpose(totvatmIO)
        s = nf90_put_var(nc, idWs,  copy_in_2d) 
       call handle_err1(s,counter,fileNetCDF)

        copy_in_2d =transpose(totempIO)
        s = nf90_put_var(nc, idE,    copy_in_2d) 
       call handle_err1(s,counter,fileNetCDF)

        copy_in_2d =transpose(totqsrIO)
        s = nf90_put_var(nc, idR,    copy_in_2d) 
       call handle_err1(s,counter,fileNetCDF)



        s= nf90_close(nc)

        deallocate(copy_in_2d)
        deallocate(copy_in)


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

       SUBROUTINE switch_index_double(in_matrix,out_matrix,x,y,z)

       IMPLICIT NONE
       integer,intent(in) :: x,y,z
       double precision,dimension(z,y,x),intent(in)    :: in_matrix
       double precision,dimension(x,y,z),intent(inout) :: out_matrix
       integer :: i,j,k

       DO i=1,x
        DO j=1,y
          DO k=1,z
            out_matrix(i,j,k) = in_matrix(k,j,i)
          ENDDO
        ENDDO
      ENDDO   

      END SUBROUTINE switch_index_double
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE switch_index_real(in_matrix,out_matrix,x,y,z)

       IMPLICIT NONE
       integer,intent(in) :: x,y,z
       real,dimension(z,y,x),intent(in)    :: in_matrix
       real,dimension(x,y,z),intent(inout) :: out_matrix
       integer :: i,j,k

       DO i=1,x
        DO j=1,y
          DO k=1,z
            out_matrix(i,j,k) = in_matrix(k,j,i)
          ENDDO
        ENDDO
      ENDDO   

      END SUBROUTINE switch_index_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      SUBROUTINE switch_index_rout(in_matrix,out_matrix,x,y,z)

       IMPLICIT NONE
       integer,intent(in) :: x,y,z
       double precision,dimension(z,y,x),intent(in)    :: in_matrix
       real,dimension(x,y,z),intent(inout) :: out_matrix
       integer :: i,j,k


       DO i=1,x
        DO j=1,y
          DO k=1,z
            out_matrix(i,j,k) = REAL(in_matrix(k,j,i),4)
          ENDDO
        ENDDO
      ENDDO   

      END SUBROUTINE switch_index_rout
