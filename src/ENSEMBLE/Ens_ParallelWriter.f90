! developed by Simone Spada (sspada@ogs.it) at OGS

module Ens_ParallelWriter
    
    use mpi
    use modul_param, &
        only: jpiglo, jpjglo, jpi, jpj, jpk
    use myalloc, &
        only: lwp, &
            nldj, nlej, nldi, nlei, &
            mycomm, mysize, myrank, &
            deflate_rst, deflate_level_rst, &
            domdec, tmask
    use Ens_Mem, &
        only: EnsDebug
        
    implicit none
    
    double precision, parameter :: PW_miss_val =1.d20
    
    character(len=1000) :: PW_filename
    character(len=100) :: PW_varname
    character(len=17) :: PW_datestring
    double precision, dimension(:,:,:), allocatable :: PW_tracer, PW_tracer_piece
    double precision, dimension(:), allocatable :: PW_tracer_puzzle
    integer, dimension(:), allocatable :: PW_displacement, PW_count, PW_displacement_2D, PW_count_2D, PW_nldi,PW_nlei,PW_nldj,PW_nlej 
    integer :: PW_n_vertical, PW_writing_rank
    double precision, dimension(:,:,:), allocatable :: nc_tracer
        
contains

subroutine PW_Init

    integer indexi, right_rank, top_rank
    
    allocate(PW_tracer(jpk, jpjglo, jpiglo))
    PW_tracer=Huge(PW_tracer(1, 1, 1))
    
    allocate(PW_nldi(mysize))
    allocate(PW_nlei(mysize))
    allocate(PW_nldj(mysize))
    allocate(PW_nlej(mysize))
    right_rank=maxval(domdec(:,2))
    top_rank = maxval(domdec(:,3))
    do indexi=1,mysize
        if (domdec(indexi, 2)==0) then
            PW_nldi(indexi)=1
        else
            PW_nldi(indexi)=2
        end if
        if (domdec(indexi, 3)==0) then
            PW_nldj(indexi)=1
        else
            PW_nldj(indexi)=2
        end if
        if (domdec(indexi, 2)==right_rank) then
            PW_nlei(indexi)=domdec(indexi,4)
        else
            PW_nlei(indexi)=domdec(indexi,4)-1
        end if
        if (domdec(indexi, 3)==top_rank) then
            PW_nlej(indexi)=domdec(indexi,5)
        else
            PW_nlej(indexi)=domdec(indexi,5)-1
        end if
    end do
    
    allocate(PW_count(mysize))
    PW_count=Huge(PW_count(1))
    
    allocate(PW_displacement(mysize+1))
    PW_displacement=Huge(PW_displacement(1))
    
    allocate(PW_count_2D(mysize))
    PW_count_2D(:) = (PW_nlej - PW_nldj +1)*(PW_nlei - PW_nldi +1)
    
    allocate(PW_displacement_2D(mysize+1))
    PW_displacement_2D(1)=0
    do indexi=1,mysize
        PW_displacement_2D(indexi+1)=PW_displacement_2D(indexi) + PW_count_2D(indexi)
    end do
    
    allocate(PW_tracer_puzzle(PW_displacement_2D(mysize+1)*jpk))
    PW_tracer_puzzle=Huge(PW_tracer_puzzle(1))
    
    PW_writing_rank=0
    
    
end subroutine

subroutine PW_Finalize
    
    deallocate(PW_tracer)    
    deallocate(PW_count)
    deallocate(PW_displacement)
    deallocate(PW_count_2D)
    deallocate(PW_displacement_2D)
    deallocate(PW_tracer_puzzle)
    deallocate(PW_nldi)
    deallocate(PW_nlei)
    deallocate(PW_nldj)
    deallocate(PW_nlej)
    
end subroutine

subroutine PW_prepare_writing(prefix, datestring, varname, tracer, n_vertical)

    character(len=17), intent(in) :: datestring
    character(len=*), intent(in) :: prefix, varname
    integer, intent(in) :: n_vertical
!     double precision, dimension(n_vertical,jpj,jpi), intent(in) :: tracer
    double precision, dimension(n_vertical, nlej - nldj +1, nlei - nldi +1), intent(in) :: tracer
    
    integer ierr, indexi, indexj, indexk
    
    allocate( PW_tracer_piece(n_vertical,nlej - nldj +1, nlei - nldi +1))
    
    !PW_tracer_piece(:,:,:)=tracer(:,nldj:nlej,nldi:nlei)*tmask(1:n_vertical,nldj:nlej,nldi:nlei) + (1-tmask(1:n_vertical,nldj:nlej,nldi:nlei))*PW_miss_val
    
    do indexi=1,nlei - nldi +1
        do indexj=1,nlej - nldj +1
            do indexk=1,n_vertical
                if (tmask(indexk,nldj+indexj-1,nldi-1+ indexi)==0) then
                    PW_tracer_piece(indexk:n_vertical,indexj,indexi)=PW_miss_val
                    exit
                end if
!                 PW_tracer_piece(indexk,indexj,indexi)=tracer(indek,nldj+indexj-1,nldi-1+ indexi)
                PW_tracer_piece(indexk,indexj,indexi)=tracer(indexk,indexj,indexi)
            end do
        end do
    end do
    
    if (myrank==PW_writing_rank) then
    
        PW_filename=trim(prefix)//'.'//trim(datestring)//'.'//trim(varname)//'.nc'
        PW_varname=trim(varname)
        PW_n_vertical=n_vertical
        PW_datestring=datestring
        
    end if
    
    PW_count=PW_count_2D*n_vertical
    PW_displacement=PW_displacement_2D*n_vertical
    CALL MPI_GATHERV(PW_tracer_piece, n_vertical*(nlej - nldj +1)*(nlei - nldi +1), MPI_DOUBLE_PRECISION, &
        PW_tracer_puzzle(1:PW_displacement(mysize+1)), PW_count, PW_displacement(1:mysize), MPI_DOUBLE_PRECISION, PW_writing_rank, mycomm, ierr)
    
    deallocate( PW_tracer_piece)
    
    if ((myrank==PW_writing_rank).and.(EnsDebug>0)) write(*,*) trim(PW_filename), " prepared for writing."
    
    PW_writing_rank=PW_writing_rank+1
    if (PW_writing_rank==mysize) call PW_write_all
    
end subroutine

subroutine PW_write_all

    integer indexi
    
    if (myrank>=PW_writing_rank) then
        PW_writing_rank=0
        return
    end if
    
    PW_count=PW_count_2D*PW_n_vertical
    PW_displacement=PW_displacement_2D*PW_n_vertical
    
    PW_tracer=PW_miss_val
    
    do indexi=1,mysize
        PW_tracer(1:PW_n_vertical, domdec(indexi,7) + PW_nldj(indexi) -1 : domdec(indexi,7) + PW_nlej(indexi) -1, domdec(indexi,6) + PW_nldi(indexi) -1 : domdec(indexi,6) + PW_nlei(indexi) -1) = &
            reshape(PW_tracer_puzzle(PW_displacement(indexi)+1:PW_count(indexi)+PW_displacement(indexi)), (/PW_n_vertical,PW_nlej(indexi) - PW_nldj(indexi) +1, PW_nlei(indexi) - PW_nldi(indexi) +1/))
    end do
    
    call PW_write(trim(PW_filename), trim(PW_varname), PW_datestring, deflate_rst, deflate_level_rst, PW_tracer, PW_n_vertical)
    
    write(*,*) trim(PW_filename), " written."
    
    PW_writing_rank=0
    
end subroutine

subroutine PW_write(fileNetCDF, VAR, TimeString, deflate, deflate_level, tracer, n_vertical)
        USE netcdf
!         USE myalloc,
!             only: totglamt, totgphit, gdept

        IMPLICIT NONE
        CHARACTER(len=*),intent(in) :: fileNetCDF
!         double precision, intent(in) :: julian
        CHARACTER(len=*),intent(in) ::  VAR
        CHARACTER(LEN=17),intent(in) :: TimeString
        integer, intent(in) :: deflate, deflate_level, n_vertical
        double precision, dimension(n_vertical,jpjglo,jpiglo), intent(in) :: tracer

        ! local
        integer :: istart, iend
        integer :: s, nc, counter
        integer :: timid, depid, yid, xid, xaid, yaid, zaid
        integer :: idB, idN, idLon, idLat, idLev, idTim
        integer shuffle
        integer indexi

        shuffle       = 0

        ! Just to try without 'or'
        ! s = nf90_create(fileNetCDF, or(nf90_clobber,NF90_HDF5), nc)
        s = nf90_create(fileNetCDF, NF90_HDF5, nc)

        s = nf90_put_att(nc, nf90_global, 'TimeString'     , TimeString)
        ! *********** DIMENSIONS ****************
        s= nf90_def_dim(nc,'x'   , jpiglo,  xid)
        s= nf90_def_dim(nc,'y'   , jpjglo,  yid)
        s= nf90_def_dim(nc,'z'   , n_vertical   ,depid)
        s= nf90_def_dim(nc,'time', 1     ,timid)

!         s = nf90_def_var(nc,'nav_lon', nf90_double,  (/xid,yid/), idLon)
!         s = nf90_def_var(nc,'nav_lat', nf90_double,  (/xid,yid/), idLat)
!         s = nf90_def_var(nc,'nav_lev', nf90_double,  (/depid/)  , idLev)
        !s = nf90_def_var(nc,'time'   , nf90_double,  (/timid/)  , idTim)
        s = nf90_def_var(nc,'TRN'//VAR, nf90_double, (/xid,yid,depid,timid/), idN)
        s = nf90_def_var_deflate(nc, idN, shuffle, deflate, deflate_level)
        call handle_err1(s,counter,fileNetCDF)
        !s= nf90_put_att(nc,idTim ,'Units', 'seconds since 1582-10-15 00:00:00')

        s = nf90_put_att(nc,idN   , 'missing_value',PW_miss_val)
        s =nf90_enddef(nc)
!         s = nf90_put_var(nc, idLon,  TRANSPOSE(totglamt))
!         call handle_err1(s,counter,fileNetCDF)
!         s = nf90_put_var(nc, idLat,  TRANSPOSE(totgphit))
!         call handle_err1(s,counter,fileNetCDF)


        !s = nf90_put_var(nc, idLev,     gdept(1:n_vertical)) !! wrong line if gdept is 3d
        !call handle_err1(s,counter,fileNetCDF)
        
        !s = nf90_put_var(nc, idN,  reshape(tracer, (/jpiglo,jpjglo,n_vertical /), ORDER= (/3,2,1/) ))
        allocate(nc_tracer(jpiglo,jpjglo,n_vertical))
        do indexi=1,jpjglo
            nc_tracer(:,indexi,:)=transpose(tracer(:,indexi,:))
        end do
        s = nf90_put_var(nc, idN, nc_tracer )
        call handle_err1(s,counter,fileNetCDF)
        s =nf90_close(nc)
        
        deallocate(nc_tracer)
        
end subroutine 

end module
