
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

       SUBROUTINE OPEN_AVE_edeseu(fileNetCDF, NLEV, VAR1, VAR2, VAR3, nc)
       USE netcdf

       IMPLICIT NONE

       integer, parameter :: nlt    = 33

       character*(*),intent(in) :: fileNetCDF

       character*(*),intent(in) :: VAR1, VAR2, VAR3

       integer,intent(in)  :: NLEV

       integer,intent(out) :: nc

       integer :: istart,iend

       integer :: ji,jj, NLEVp1

       integer :: s, counter

       integer :: time_id, depth_id, nid

       integer :: idVAR1, idVAR2, idVAR3

        s = nf90_create(fileNetCDF, or(nf90_clobber,NF90_HDF5), nc)
        ! *********** GLOBAL ********************
        s = nf90_put_att(nc, nf90_global, 'Convenctions'   ,'COARDS')

        ! *********** DIMENSIONS ****************
        s= nf90_def_dim(nc,'nl'  , nlt,   nid)

        NLEVp1 = NLEV +1

        s= nf90_def_dim(nc,'depth', NLEVp1, depth_id)

        ! ********** VARIABLES *****************

        s = nf90_def_var(nc,trim(VAR1) , nf90_double, (/depth_id, nid/), idVAR1)

        s = nf90_def_var(nc,trim(VAR2) , nf90_double, (/depth_id, nid/), idVAR2) 

        s = nf90_def_var(nc,trim(VAR3) , nf90_double, (/depth_id, nid/), idVAR3)


        s = nf90_put_att(nc,idVAR1, 'long_name'    ,trim(VAR1))

        s = nf90_put_att(nc,idVAR1, 'missing_value' ,1.e+20)

        s = nf90_put_att(nc,idVAR2, 'long_name'    ,trim(VAR2))

        s = nf90_put_att(nc,idVAR2, 'missing_value' ,1.e+20)

        s = nf90_put_att(nc,idVAR3, 'long_name'    ,trim(VAR3))

        s = nf90_put_att(nc,idVAR3, 'missing_value' ,1.e+20)

        s =nf90_enddef(nc)

       END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       SUBROUTINE WRITE_AVE_edeseu(fileNetCDF, nc, NLEV, VAR1, VAR2, VAR3, M1_0, M2_0, M3_0, M1, M2, M3)
       USE netcdf
       IMPLICIT NONE
       integer, parameter :: nlt    = 33

       character*(*),intent(in) :: fileNetCDF

       integer, intent(in)      :: nc, NLEV

       character*(*),intent(in) :: VAR1, VAR2, VAR3

       double precision ,intent(in),dimension(NLEV, nlt) :: M1, M2, M3

       double precision ,intent(in),dimension(nlt) :: M1_0, M2_0, M3_0

       double precision, allocatable :: M1p1(:,:), M2p1(:,:), M3p1(:,:)

       integer :: NLEVp1,s,idVAR, counter

       integer :: thestart(2), thecount(2)

        NLEVp1 = NLEV + 1

        allocate(M1p1(NLEVp1,nlt))
        allocate(M2p1(NLEVp1,nlt))
        allocate(M3p1(NLEVp1,nlt))

        M1p1(1,:) = M1_0
        M2p1(1,:) = M2_0
        M3p1(1,:) = M3_0


        M1p1(2:NLEVp1,:) = M1
        M2p1(2:NLEVp1,:) = M2
        M3p1(2:NLEVp1,:) = M3

        thestart = (/1, 1 /)
        thecount = (/NLEVp1, nlt/)

        counter = 0

        s = nf90_inq_varid (nc, trim(VAR1), idVAR) ; call handle_err1(s,fileNetCDF)
        s = nf90_put_var(nc, idVAR ,  M1p1, thestart, thecount) 
         call handle_err1(s,counter,fileNetCDF)

        s = nf90_inq_varid (nc, trim(VAR2), idVAR) ; call handle_err1(s,fileNetCDF)
        s = nf90_put_var(nc, idVAR ,  M2p1, thestart, thecount) 
         call handle_err1(s,counter,fileNetCDF)

        s = nf90_inq_varid (nc, trim(VAR3), idVAR) ; call handle_err1(s,fileNetCDF)
        s = nf90_put_var(nc, idVAR ,  M3p1, thestart, thecount) 
         call handle_err1(s,counter,fileNetCDF)

       END SUBROUTINE

       SUBROUTINE CLOSE_AVE_edeseu(fileNetCDF,nc)
       USE netcdf
       IMPLICIT NONE

       CHARACTER*(*), INTENT(IN)    :: fileNetCDF
       integer,       INTENT(IN)    :: nc
       ! Local variable
       integer                      :: s

        s =nf90_close(nc); call handle_err1(s,fileNetCDF)

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
