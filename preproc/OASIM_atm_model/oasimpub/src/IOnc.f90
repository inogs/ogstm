
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


! ************************************************************

      SUBROUTINE readnc_global_float_2d(fileNetCDF,varname, M)
      USE netcdf

      implicit none
      
      character,intent(in)    :: fileNetCDF*(*) ,varname*(*)
      REAL     ,intent(inout) :: M(360,180)

      integer ncid, stat, VARid,i,j
      integer counter
      integer thecount(2), start(2)
      REAL    :: INDATA(360,180)

      counter = 0
      
      start    = (/1,       1/)
      thecount = (/360,   180/)

      stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)  
      call handle_err1(stat, counter,FileNetCDF)

      stat = nf90_inq_varid (ncid, varname, VARid)
      call handle_err2(stat, fileNetCDF,varname)
      call handle_err1(stat, counter,FileNetCDF)

      stat = nf90_get_var (ncid,VARid,INDATA,start, thecount)
      call handle_err2(stat, fileNetCDF,varname)        
      call handle_err1(stat, counter,FileNetCDF)

      stat = nf90_close(ncid)                           
      call handle_err1(stat, counter,FileNetCDF)

      M(:,:) = INDATA(:,:)

      END SUBROUTINE readnc_global_float_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      INTERFACE WRITE_AVE_2D

!      SUBROUTINE WRITE_AVE_2D_1v(fileNetCDF, VAR1, M1)
!         integer, parameter :: jpiglo = 360 
!         integer, parameter :: jpjglo = 180
!         character*(*),intent(in) :: fileNetCDF
!         character(LEN=20),intent(in) :: VAR1
!         real,intent(in),dimension(jpiglo, jpjglo) :: M1
!      END SUBROUTINE WRITE_AVE_2D_1v

!      SUBROUTINE WRITE_AVE_2D_2v_l(fileNetCDF, VAR1, VAR2, M1, M2)
!         integer, parameter :: jpiglo = 360 
!         integer, parameter :: jpjglo = 180
!         integer, parameter :: nlt    = 33
!         character*(*),intent(in) :: fileNetCDF
!         character(LEN=20),intent(in) :: VAR1,VAR2
!         real,intent(in),dimension(jpiglo, jpjglo, nlt) :: M1, M2
!      END SUBROUTINE WRITE_AVE_2D_2v_l

!      SUBROUTINE WRITE_AVE_2D_2v_l_h(fileNetCDF, VAR1, VAR2, M1, M2)
!         integer, parameter :: jpiglo = 360 
!         integer, parameter :: jpjglo = 180
!         integer, parameter :: nlt    = 33
!         integer, parameter :: ntime  = 12
!         character*(*),intent(in) :: fileNetCDF
!         character(LEN=20),intent(in) :: VAR1,VAR2
!         real,intent(in),dimension(jpiglo, jpjglo, nlt, ntime) :: M1, M2
!      END SUBROUTINE WRITE_AVE_2D_2v_l_h

!      END INTERFACE

       SUBROUTINE WRITE_AVE_2D_1v(fileNetCDF, VAR1, M1)
       USE netcdf
       IMPLICIT NONE
       integer, parameter :: jpiglo = 360 
       integer, parameter :: jpjglo = 180

       character*(*),intent(in) :: fileNetCDF
       character*(*),intent(in) :: VAR1
       real,intent(in),dimension(jpiglo, jpjglo) :: M1

       integer :: istart,iend

       integer :: ji,jj
       integer :: s, nc, counter
       integer :: timid, yid, xid, nid
       integer :: idvartime, idphit, idlamt, idVAR1

       real :: lat_actual_range(2), lon_actual_range(2)
       REAL :: totglamt(jpiglo,jpjglo),  totgphit(jpiglo,jpjglo)
         lon_actual_range=(/-180.0 , 180.0   /)
         lat_actual_range=(/-90.0  , 90.0    /)


        do jj=1,jpjglo
            do ji=1,jpiglo
                totglamt(ji,jj) = -180.0 + REAL(ji) -1.0
                totgphit(ji,jj) = -90.0  + REAL(jj) -1.0
            enddo
        enddo

        s = nf90_create(fileNetCDF, or(nf90_clobber,NF90_HDF5), nc)
        ! *********** GLOBAL ********************
        s = nf90_put_att(nc, nf90_global, 'Convenctions'   ,'COARDS')

        ! *********** DIMENSIONS ****************
        s= nf90_def_dim(nc,'lon'           , jpiglo,  xid)
        s= nf90_def_dim(nc,'lat'           , jpjglo,  yid)

        ! ********** VARIABLES *****************
        s = nf90_def_var(nc,'lat'   ,       nf90_float, (/yid/), idphit)
        s = nf90_def_var(nc,'lon'   ,       nf90_float, (/xid/), idlamt)

        s = nf90_def_var(nc,trim(VAR1) ,    nf90_float, (/xid,yid/),  idVAR1)


        s = nf90_put_att(nc,idphit, 'units'        ,'degrees_north')
        s = nf90_put_att(nc,idphit, 'long_name'    ,'Latitude')
        s = nf90_put_att(nc,idphit,'actual_range' ,lat_actual_range)

        s = nf90_put_att(nc,idlamt, 'units'        ,'degrees_east')
        s = nf90_put_att(nc,idlamt, 'long_name'    ,'Longitude')
        s = nf90_put_att(nc,idlamt,'actual_range' ,lon_actual_range)

        s = nf90_put_att(nc,idVAR1, 'long_name'    ,trim(VAR1))
        s = nf90_put_att(nc,idVAR1, 'missing_value' ,1.e+20)

        s =nf90_enddef(nc)

        counter=0
        s = nf90_put_var(nc, idlamt,  REAL(totglamt(:,jpjglo),4) )
        call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idphit,  REAL(totgphit(jpiglo,:),4) )
        call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idVAR1  ,  M1) 
        call handle_err1(s,counter,fileNetCDF)

        s =nf90_close(nc)

       END SUBROUTINE 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       SUBROUTINE WRITE_AVE_2D_2v_l(fileNetCDF, VAR1, VAR2, M1, M2)
       USE netcdf
       IMPLICIT NONE
       integer, parameter :: jpiglo = 360 
       integer, parameter :: jpjglo = 180
       integer, parameter :: nlt    = 33

       character*(*),intent(in) :: fileNetCDF
       character*(*),intent(in) :: VAR1, VAR2
       real,intent(in),dimension(jpiglo, jpjglo, nlt) :: M1, M2

       integer :: istart,iend

       integer :: ji,jj
       integer :: s, nc, counter
       integer :: yid, xid, nid
       integer :: idvartime, idphit, idlamt, idVAR1, idVAR2

       real :: lat_actual_range(2), lon_actual_range(2)
       REAL :: totglamt(jpiglo,jpjglo),  totgphit(jpiglo,jpjglo)
         lon_actual_range=(/-180.0 , 180.0   /)
         lat_actual_range=(/-90.0  , 90.0    /)


        do jj=1,jpjglo
            do ji=1,jpiglo
                totglamt(ji,jj) = -180.0 + REAL(ji) -1.0
                totgphit(ji,jj) = -90.0  + REAL(jj) -1.0
            enddo
        enddo

        s = nf90_create(fileNetCDF, or(nf90_clobber,NF90_HDF5), nc)
        ! *********** GLOBAL ********************
        s = nf90_put_att(nc, nf90_global, 'Convenctions'   ,'COARDS')

        ! *********** DIMENSIONS ****************
        s= nf90_def_dim(nc,'lon'           , jpiglo,  xid)
        s= nf90_def_dim(nc,'lat'           , jpjglo,  yid)
        s= nf90_def_dim(nc,'nl'            ,    nlt,  nid)

        ! ********** VARIABLES *****************
        !s = nf90_def_var(nc,'time',         nf90_double,(/timid/),
        !idvartime)
        s = nf90_def_var(nc,'lat'   ,       nf90_float, (/yid/), idphit)
        s = nf90_def_var(nc,'lon'   ,       nf90_float, (/xid/), idlamt)

        s = nf90_def_var(nc,trim(VAR1) ,        nf90_float, (/xid,yid,nid/),  idVAR1)

        s = nf90_def_var(nc,trim(VAR2) ,        nf90_float, (/xid,yid,nid/),  idVAR2)

        s = nf90_put_att(nc,idphit, 'units'        ,'degrees_north')
        s = nf90_put_att(nc,idphit, 'long_name'    ,'Latitude')
        s = nf90_put_att(nc,idphit,'actual_range' ,lat_actual_range)

        s = nf90_put_att(nc,idlamt, 'units'        ,'degrees_east')
        s = nf90_put_att(nc,idlamt, 'long_name'    ,'Longitude')
        s = nf90_put_att(nc,idlamt,'actual_range' ,lon_actual_range)

        s = nf90_put_att(nc,idVAR1, 'long_name'    ,trim(VAR1))
        s = nf90_put_att(nc,idVAR1, 'missing_value' ,1.e+20)

        s = nf90_put_att(nc,idVAR2, 'long_name'    ,trim(VAR2))
        s = nf90_put_att(nc,idVAR2, 'missing_value' ,1.e+20)

        s =nf90_enddef(nc)

        counter=0
        s = nf90_put_var(nc, idlamt,  REAL(totglamt(:,jpjglo),4) )
        call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idphit,  REAL(totgphit(jpiglo,:),4) )
        call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idVAR1  ,  M1) 
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idVAR2 ,  M2) 
       call handle_err1(s,counter,fileNetCDF)


        s =nf90_close(nc)

       END SUBROUTINE 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       SUBROUTINE WRITE_AVE_2D_2v_l_h(fileNetCDF, VAR1, VAR2, M1, M2)
!      SUBROUTINE WRITE_AVE_2D(fileNetCDF, VAR1, VAR2, datefrom, dateTo, M1, M2)
       USE netcdf
       IMPLICIT NONE
       integer, parameter :: jpiglo = 360 
       integer, parameter :: jpjglo = 180
       integer, parameter :: nlt    = 33
       integer, parameter :: ntime  = 12

       character*(*),intent(in) :: fileNetCDF
       character*(*),intent(in) :: VAR1, VAR2
       real,intent(in),dimension(jpiglo, jpjglo, nlt, ntime) :: M1, M2

       integer :: istart,iend

       integer :: ji,jj
       integer :: s, nc, counter
       integer :: timid, yid, xid, nid
       integer :: idvartime, idphit, idlamt, idVAR1, idVAR2

       real :: lat_actual_range(2), lon_actual_range(2)
       REAL :: totglamt(jpiglo,jpjglo),  totgphit(jpiglo,jpjglo)
         lon_actual_range=(/-180.0 , 180.0   /)
         lat_actual_range=(/-90.0  , 90.0    /)


        do jj=1,jpjglo
            do ji=1,jpiglo
                totglamt(ji,jj) = -180.0 + REAL(ji) -1.0
                totgphit(ji,jj) = -90.0  + REAL(jj) -1.0
            enddo
        enddo

        s = nf90_create(fileNetCDF, or(nf90_clobber,NF90_HDF5), nc)
        ! *********** GLOBAL ********************
        s = nf90_put_att(nc, nf90_global, 'Convenctions'   ,'COARDS')
!       s = nf90_put_att(nc, nf90_global, 'DateStart'     , datefrom)
!       s = nf90_put_att(nc, nf90_global, 'Date__End'     ,   dateTo)

        ! *********** DIMENSIONS ****************
        s= nf90_def_dim(nc,'lon'           , jpiglo,  xid)
        s= nf90_def_dim(nc,'lat'           , jpjglo,  yid)
        s= nf90_def_dim(nc,'nl'            ,    nlt,  nid)
        s= nf90_def_dim(nc,'time'          ,  ntime,timid)

        ! ********** VARIABLES *****************
        !s = nf90_def_var(nc,'time',         nf90_double,(/timid/),       idvartime)
        s = nf90_def_var(nc,'lat'   ,       nf90_float, (/yid/),            idphit)
        s = nf90_def_var(nc,'lon'   ,       nf90_float, (/xid/),            idlamt)

        s = nf90_def_var(nc,trim(VAR1) ,        nf90_float, (/xid,yid,nid,timid/),  idVAR1)

        s = nf90_def_var(nc,trim(VAR2) ,        nf90_float, (/xid,yid,nid,timid/),  idVAR2)

        s = nf90_put_att(nc,idphit, 'units'        ,'degrees_north')
        s = nf90_put_att(nc,idphit, 'long_name'    ,'Latitude')
        s = nf90_put_att(nc,idphit,'actual_range' ,lat_actual_range)

        s = nf90_put_att(nc,idlamt, 'units'        ,'degrees_east')
        s = nf90_put_att(nc,idlamt, 'long_name'    ,'Longitude')
        s = nf90_put_att(nc,idlamt,'actual_range' ,lon_actual_range)

        s = nf90_put_att(nc,idVAR1, 'long_name'    ,trim(VAR1))
        s = nf90_put_att(nc,idVAR1, 'missing_value' ,1.e+20)

        s = nf90_put_att(nc,idVAR2, 'long_name'    ,trim(VAR2))
        s = nf90_put_att(nc,idVAR2, 'missing_value' ,1.e+20)

        s =nf90_enddef(nc)

        counter=0
        s = nf90_put_var(nc, idlamt,  REAL(totglamt(:,jpjglo),4) )
        call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idphit,  REAL(totgphit(jpiglo,:),4) )
        call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idVAR1  ,  M1) 
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idVAR2 ,  M2) 
       call handle_err1(s,counter,fileNetCDF)


        s =nf90_close(nc)

       END SUBROUTINE


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
