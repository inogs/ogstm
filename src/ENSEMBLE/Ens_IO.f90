! developed by Simone Spada (sspada@ogs.it) at OGS

module Ens_IO
    use modul_param, &
        only: jpiglo, jpjglo, jpi, jpj, jpk, jptra, jptra_dia, jptra_dia_2d
    use myalloc, &
        only: jptra_high, jptra_dia_high, jptra_dia2d_high, jptra_phys, jptra_phys_2d, &
            traIO, traIO_HIGH, tra_DIA_IO, tra_DIA_IO_HIGH, tra_DIA_2d_IO, tra_DIA_2d_IO_HIGH, &
            tra_PHYS_IO, tra_PHYS_IO_HIGH, tra_PHYS_2d_IO, tra_PHYS_2d_IO_HIGH, &
            freq_ave_phys, snIO, tnIO, vatmIO, empIO, qsrIO, unIO, vnIO, wnIO, avtIO, e3tIO
    use DIA_mem, &
        only: Fsize, diaflx
    use MPI_GATHER_INFO, &
            only: WRITING_RANK_WR
            
    implicit none
    double precision, dimension(:,:,:), allocatable :: Ens_copy_in, Ens_copy_in_real
    integer :: n_traIO, n_traIO_HIGH, n_tra_DIA_IO, n_tra_DIA_IO_HIGH, n_tra_DIA_2d_IO, n_tra_DIA_2d_IO_HIGH, &
        n_tra_PHYS_IO, n_tra_PHYS_IO_HIGH, n_tra_PHYS_2d_IO, n_tra_PHYS_2d_IO_HIGH, n_diaflx
    integer :: win_traIO, win_traIO_HIGH, win_tra_DIA_IO, win_tra_DIA_IO_HIGH, win_tra_DIA_2d_IO, win_tra_DIA_2d_IO_HIGH, &
        win_tra_PHYS_IO, win_tra_PHYS_IO_HIGH, win_tra_PHYS_2d_IO, win_tra_PHYS_2d_IO_HIGH, win_diaflx
    double precision, dimension(:,:), contiguous, pointer :: gl_traIO, gl_traIO_HIGH, gl_tra_DIA_IO, gl_tra_DIA_IO_HIGH, gl_tra_DIA_2d_IO, gl_tra_DIA_2d_IO_HIGH, &
        gl_tra_PHYS_IO, gl_tra_PHYS_IO_HIGH, gl_tra_PHYS_2d_IO, gl_tra_PHYS_2d_IO_HIGH, gl_diaflx
        
contains

    Subroutine Ens_Init_IO
        use Ens_Mem, &
            only: EnsSize, EnsAveDouble, & 
                Ens_shared_alloc
        
        double precision, dimension(:), POINTER, contiguous :: member_pointer
        double precision, dimension(:,:), POINTER, contiguous :: global_pointer
        
        if (WRITING_RANK_WR) then
            allocate(Ens_copy_in(jpiglo, jpjglo, jpk))
            if (.not.EnsAveDouble) allocate(Ens_copy_in_real(jpiglo, jpjglo, jpk))
        end if
        
        
        n_traIO=jpk*jpj*jpi*jptra
        call Ens_shared_alloc(n_traIO, member_pointer, global_pointer, win_traIO)
        traIO(1:jpk,1:jpj,1:jpi,1:jptra)=>member_pointer
        gl_traIO(1:n_traIO, 0:EnsSize-1)=>global_pointer
        traIO  = huge(traIO(1,1,1,1)) 
        
        n_traIO_HIGH=jpk*jpj*jpi*jptra_HIGH
        call Ens_shared_alloc(n_traIO_HIGH, member_pointer, global_pointer, win_traIO_HIGH)
        traIO_HIGH(1:jpk,1:jpj,1:jpi,1:jptra_HIGH)=>member_pointer
        gl_traIO_HIGH(1:n_traIO_HIGH, 0:EnsSize-1)=>global_pointer
        traIO_HIGH    = huge(traIO_HIGH(1,1,1,1))
        
        
        n_tra_DIA_IO=jpk*jpj*jpi*jptra_dia
        call Ens_shared_alloc(n_tra_DIA_IO, member_pointer, global_pointer, win_tra_DIA_IO)
        tra_DIA_IO(1:jptra_dia,1:jpk,1:jpj,1:jpi)=>member_pointer
        gl_tra_DIA_IO(1:n_tra_DIA_IO, 0:EnsSize-1)=>global_pointer
        tra_DIA_IO    = huge(tra_DIA_IO(1,1,1,1))
        
        n_tra_DIA_IO_HIGH=jpk*jpj*jpi*jptra_dia_HIGH
        call Ens_shared_alloc(n_tra_DIA_IO_HIGH, member_pointer, global_pointer, win_tra_DIA_IO_HIGH)
        tra_DIA_IO_HIGH(1:jptra_dia_HIGH,1:jpk,1:jpj,1:jpi)=>member_pointer
        gl_tra_DIA_IO_HIGH(1:n_tra_DIA_IO_HIGH, 0:EnsSize-1)=>global_pointer
        tra_DIA_IO_HIGH = huge(tra_DIA_IO_HIGH(1,1,1,1))
        
        n_tra_DIA_2d_IO=jpj*jpi*jptra_dia_2d
        call Ens_shared_alloc(n_tra_DIA_2d_IO, member_pointer, global_pointer, win_tra_DIA_2d_IO)
        tra_DIA_2d_IO(1:jptra_dia_2d,1:jpj,1:jpi)=>member_pointer
        gl_tra_DIA_2d_IO(1:n_tra_DIA_2d_IO, 0:EnsSize-1)=>global_pointer
        tra_DIA_2d_IO    = huge(tra_DIA_2d_IO(1,1,1))
    
        n_tra_DIA_2d_IO_HIGH=jpj*jpi*jptra_dia2d_HIGH
        call Ens_shared_alloc(n_tra_DIA_2d_IO_HIGH, member_pointer, global_pointer, win_tra_DIA_2d_IO_HIGH)
        tra_DIA_2d_IO_HIGH(1:jptra_dia2d_HIGH,1:jpj,1:jpi)=>member_pointer
        gl_tra_DIA_2d_IO_HIGH(1:n_tra_DIA_2d_IO_HIGH, 0:EnsSize-1)=>global_pointer
        tra_DIA_2d_IO_HIGH = huge(tra_DIA_2d_IO_HIGH(1,1,1))
        
        
        if (freq_ave_phys==1) then
            n_tra_PHYS_IO_HIGH=jpk*jpj*jpi*jptra_phys
            call Ens_shared_alloc(n_tra_PHYS_IO_HIGH, member_pointer, global_pointer, win_tra_PHYS_IO_HIGH)
            tra_PHYS_IO_HIGH(1:jpk,1:jpj,1:jpi,1:jptra_phys)=>member_pointer
            gl_tra_PHYS_IO_HIGH(1:n_tra_PHYS_IO_HIGH, 0:EnsSize-1)=>global_pointer
            tra_PHYS_IO_HIGH = huge(tra_PHYS_IO_HIGH(1,1,1,1))
            
            n_tra_PHYS_2d_IO_HIGH=jpj*jpi*jptra_phys_2d
            call Ens_shared_alloc(n_tra_PHYS_2d_IO_HIGH, member_pointer, global_pointer, win_tra_PHYS_2d_IO_HIGH)
            tra_PHYS_2d_IO_HIGH(1:jpj,1:jpi,1:jptra_phys_2d)=>member_pointer
            gl_tra_PHYS_2d_IO_HIGH(1:n_tra_PHYS_2d_IO_HIGH, 0:EnsSize-1)=>global_pointer
            tra_PHYS_2d_IO_HIGH = huge(tra_PHYS_2d_IO_HIGH(1,1,1))
        end if
        
        if (freq_ave_phys==2) then
            n_tra_PHYS_IO=jpk*jpj*jpi*jptra_phys
            call Ens_shared_alloc(n_tra_PHYS_IO, member_pointer, global_pointer, win_tra_PHYS_IO)
            tra_PHYS_IO(1:jpk,1:jpj,1:jpi,1:jptra_phys)=>member_pointer
            gl_tra_PHYS_IO(1:n_tra_PHYS_IO, 0:EnsSize-1)=>global_pointer
            tra_PHYS_IO    = huge(tra_PHYS_IO(1,1,1,1))
            
            n_tra_PHYS_2d_IO=jpj*jpi*jptra_phys_2d
            call Ens_shared_alloc(n_tra_PHYS_2d_IO, member_pointer, global_pointer, win_tra_PHYS_2d_IO)
            tra_PHYS_2d_IO(1:jpj,1:jpi,1:jptra_phys_2d)=>member_pointer
            gl_tra_PHYS_2d_IO(1:n_tra_PHYS_2d_IO, 0:EnsSize-1)=>global_pointer
            tra_PHYS_2d_IO    = huge(tra_PHYS_2d_IO(1,1,1))
        end if
        
        
        if (Fsize.NE.0) then
            n_diaflx=7*Fsize*jptra
            call Ens_shared_alloc(n_diaflx, member_pointer, global_pointer, win_diaflx)
            diaflx(1:7, 1:Fsize, 1:jptra)=>member_pointer
            gl_diaflx(1:n_diaflx, 0:EnsSize-1)=>global_pointer
            diaflx    = 0.0d0
        end if
        
    end subroutine
    
    Subroutine Ens_Finalize_IO
        use mpi
        
        use Ens_Mem, &
            only: EnsAveDouble
        
        integer ierror
        
        if (WRITING_RANK_WR) then
            deallocate(Ens_copy_in)          
            if (.not.EnsAveDouble) deallocate(Ens_copy_in_real)          
        end if        
        
        CALL MPI_Win_free(win_traIO, ierror)
        CALL MPI_Win_free(win_traIO_HIGH, ierror)
        
        CALL MPI_Win_free(win_tra_DIA_IO, ierror)
        CALL MPI_Win_free(win_tra_DIA_IO_HIGH, ierror)
        CALL MPI_Win_free(win_tra_DIA_2d_IO, ierror)
        CALL MPI_Win_free(win_tra_DIA_2d_IO_HIGH, ierror)
        
        if (freq_ave_phys==1) then
            CALL MPI_Win_free(win_tra_PHYS_IO_HIGH, ierror)
            CALL MPI_Win_free(win_tra_PHYS_2d_IO_HIGH, ierror)
        end if
        if (freq_ave_phys==2) then
            CALL MPI_Win_free(win_tra_PHYS_IO, ierror)
            CALL MPI_Win_free(win_tra_PHYS_2d_IO, ierror)
        end if
        
        if (Fsize.NE.0) CALL MPI_Win_free(win_diaflx, ierror)
        
    end subroutine

    SUBROUTINE Ens_trcwri(prefix, TimeString, BaseIndex, Tracer)
        use Ens_Utilities, &
            only: int2str
        
        CHARACTER(LEN=17), INTENT(IN) :: TimeString
        integer, intent(in) :: BaseIndex
        CHARACTER(len=*), intent(in) :: prefix
        double precision, dimension(jpk,jpj,jpi,jptra), intent(in) :: Tracer
        
        if (BaseIndex==-1) then
            call trcwri_Wrapped(prefix, TimeString, Tracer)
        else
            call trcwri_Wrapped(prefix//int2str(BaseIndex,3) , TimeString, Tracer)
        end if
    end SUBROUTINE

    SUBROUTINE trcwri_Wrapped(prefix, datestring, tracer)

        USE myalloc
        USE IO_mem
        USE calendar
        USE TIME_MANAGER
        use mpi
        USE ogstm_mpi_module

        USE MPI_GATHER_INFO

        USE MATRIX_VARS
        USE NODES_MODULE
        USE DTYPE_PROCS_STRING_MODULE

        IMPLICIT NONE
        CHARACTER(LEN=17), INTENT(IN) :: datestring
        CHARACTER(len=*), intent(in) :: prefix
        double precision, dimension(jpk,jpj,jpi,jptra), intent(in) :: Tracer

!----------------------------------------------------------------------
! local declarations
! ==================
        double precision ::  Miss_val =1.e20
        INTEGER jk,jj,ji,jn
        double precision julian

        CHARACTER(LEN=3) varname

        INTEGER idrank, ierr, istart, jstart, iPe, iPd, jPe, jPd, status(MPI_STATUS_SIZE)
        INTEGER irange, jrange
        INTEGER totistart, totiend, relistart, reliend
        INTEGER totjstart, totjend, reljstart, reljend
        INTEGER ind1, i_contribution, j_contribution
        CHARACTER(LEN=20)  var_to_store
        INTEGER :: COUNTER_VAR_TRCWRI, n_dumping_cycles, jv, ivar, writing_rank
        
        if (MyRank==0) write(*,*) 'Writing ',prefix

        !filename = 'RST.20111231-15:30:00.N1p.nc'
        julian=datestring2sec(datestring)

        trcwriparttime = MPI_WTIME() ! cronometer-start

        bufftrn = Miss_val
        if (WRITING_RANK_WR) tottrn = Miss_val


        n_dumping_cycles = matrix_state_2_row

        COUNTER_VAR_TRCWRI = 1

        RESTARTS_LOOP: DO jv = 1, n_dumping_cycles
                                
                DO ivar = 1 , nodes
        
                        writing_rank = writing_procs(ivar)
                        IF (COUNTER_VAR_TRCWRI > JPTRA)then
                                EXIT
                        ELSE

                                do ji =1 , jpi
                                        i_contribution= jpk*jpj * (ji - 1 )
                                        do jj =1 , jpj
                                                j_contribution=jpk*(jj-1)
                                                do jk =1 , jpk
                                                        ind1 = jk + j_contribution + i_contribution
                                                        if (tmask(jk,jj,ji).eq.1) then
                                                                bufftrn(ind1)= Tracer(jk,jj,ji,COUNTER_VAR_TRCWRI)
                                                        endif
                                                enddo
                                        enddo
                                enddo
                                COUNTER_VAR_TRCWRI = COUNTER_VAR_TRCWRI + 1

                                CALL MPI_GATHERV(bufftrn, sendcount, MPI_DOUBLE_PRECISION, bufftrn_TOT,jprcv_count, jpdispl_count,MPI_DOUBLE_PRECISION, writing_rank,mycomm, IERR)

                        END IF

                END DO

                if(WRITING_RANK_WR) then
                
                        var_to_store =matrix_state_2(jv,ind_col)%var_name
                        IF (var_to_store == "novars_input")then
                                EXIT
                        ELSE
                                DO idrank = 0,mysize-1
                                        ! ******* WRITING RANK sets
                                        ! indexes of tot matrix where to
                                        ! place buffers of idrank
                                        irange    = iPe_a(idrank+1) -iPd_a(idrank+1) + 1
                                        jrange    = jPe_a(idrank+1) -jPd_a(idrank+1) + 1
                                        totistart = istart_a(idrank+1) +iPd_a(idrank+1) - 1
                                        totiend   = totistart + irange -1
                                        totjstart = jstart_a(idrank+1) +jPd_a(idrank+1) - 1
                                        totjend   = totjstart + jrange -1
                                        relistart = 1 + iPd_a(idrank+1)- 1
                                        reliend   = relistart + irange -1
                                        reljstart = 1 + jPd_a(idrank+1)- 1
                                        reljend   = reljstart + jrange -1
                                        ! **** ASSEMBLING *** WRITING
                                        ! RANK  puts in tot matrix
                                        ! buffer received by idrank
                                        do ji =totistart,totiend
                                                i_contribution   =jpk*jpj_rec_a(idrank+1)*(ji-1-totistart+ relistart)
                                                do jj =totjstart,totjend
                                                        j_contribution =jpk*(jj-1-totjstart+ reljstart)
                                                        do jk =1, jpk
                                                                ind1 =jk + j_contribution + i_contribution
                                                                tottrn(jk,jj,ji)=bufftrn_TOT(ind1+jpdispl_count(idrank+1))
                                                        enddo
                                                enddo
                                        enddo
                                END DO

                                CALL Ens_write_restart(prefix//'.'//datestring//'.'//trim(var_to_store)//'.nc', &
                                    var_to_store, datestring, julian, deflate_rst, deflate_level_rst)

                        END IF
                END IF
        END DO RESTARTS_LOOP
    END SUBROUTINE

    SUBROUTINE Ens_write_restart(fileNetCDF,VAR,TimeString, julian,deflate, deflate_level)
        USE netcdf
        USE myalloc

        IMPLICIT NONE
        CHARACTER*(*),intent(in) :: fileNetCDF
        double precision,intent(in) :: julian
        CHARACTER(*),intent(in) ::  VAR
        CHARACTER(LEN=17),intent(in) :: TimeString
        integer, intent(in) :: deflate, deflate_level

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

        !call switch_index_double(tottrn,copy_in,jpiglo,jpjglo,jpk)
        do indexi=1, jpjglo
            Ens_copy_in(:,indexi,:)=transpose(tottrn(:,indexi,:))
        end do

        call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idN,      Ens_copy_in)
        call handle_err1(s,counter,fileNetCDF)
        s =nf90_close(nc)


    END SUBROUTINE 
    
    SUBROUTINE Ens_WRITE_AVE(fileNetCDF,VAR, datefrom, dateTo,M,deflate, deflate_level)
       USE netcdf
!use mpi     

       USE myalloc
       Use Ens_Mem, &
            only: EnsAveDouble

       IMPLICIT NONE

       CHARACTER*(*),intent(in) :: fileNetCDF
       character(LEN=20),intent(in) :: VAR
       character(LEN=17),intent(in) :: datefrom, dateTo
       double precision, dimension(jpk, jpjglo, jpiglo),intent(in) :: M
       integer, intent(in) :: deflate, deflate_level
       
       integer istart,iend
       integer s, nc, counter
       integer timid, depid, yid, xid
       integer idvartime,idgdept,idphit,idlamt,idVAR
       integer shuffle
       real lat_actual_range(2), lon_actual_range(2), depth_actual_range(2)
       integer indexi
       
!double precision tempo

!tempo=mpi_wtime()
!real, allocatable :: ciccio(:, :, :)
!write(*,*) 'comincio. t: ', mpi_wtime()-tempo
!allocate(ciccio(jpiglo, jpjglo, jpk))
!write(*,*) 'allocato'
       
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
        !s= nf90_def_dim(nc,'time'  , NF90_UNLIMITED,timid)
        s= nf90_def_dim(nc,'time'  , 1,timid)
       call handle_err1(s,counter,fileNetCDF)
       
!write(*,*) 'definisco le variabili. t: ', mpi_wtime()-tempo
        ! ********** VARIABLES *****************
        !!s = nf90_def_var(nc,'time',         nf90_double,(/timid/),       idvartime)
        !call handle_err1(s,counter,fileNetCDF)
        s = nf90_def_var(nc,'depth',        nf90_float, (/depid/),         idgdept)
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_def_var(nc,'lat'   ,       nf90_float, (/yid/),            idphit)
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_def_var(nc,'lon'   ,       nf90_float, (/xid/),            idlamt)
       call handle_err1(s,counter,fileNetCDF)
       
!write(*,*) 'definisco quella grossa. t: ', mpi_wtime()-tempo  
        if (EnsAveDouble) then
            s = nf90_def_var(nc,trim(VAR) ,        nf90_double, (/xid,yid,depid,timid/),  idVAR)
        else
            s = nf90_def_var(nc,trim(VAR) ,        nf90_float, (/xid,yid,depid,timid/),  idVAR)
        end if
        
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

!write(*,*) 'aggiungo i valori t: ', mpi_wtime()-tempo
        s = nf90_put_var(nc, idlamt,   REAL(totglamt(jpjglo,:),4) )
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idphit,   REAL(totgphit(:,jpiglo),4) )
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idgdept,  REAL(   gdept,     4) )
       call handle_err1(s,counter,fileNetCDF)
       
!write(*,*) 'prima dello scambio indici t: ', mpi_wtime()-tempo
       !call switch_index_rout(M,copy_in,jpiglo,jpjglo,jpk)
       if (EnsAveDouble) then
            do indexi=1, jpjglo
                Ens_copy_in(:,indexi,:)=transpose(M(:,indexi,:))
            end do
        else
            do indexi=1, jpjglo
                Ens_copy_in_real(:,indexi,:)=transpose(M(:,indexi,:))
            end do
        end if
            
!write(*,*) 'sum', sum(Ens_copy_in)
!write(*,*) 'sum', sum(real(Ens_copy_in,4))
!write(*,*) 'prima del casting'
!ciccio=sngl(Ens_copy_in)
!write(*,*) 'dopo il casting t: ', mpi_wtime()-tempo 
!s = nf90_put_var(nc, idVAR  ,  sngl(Ens_copy_in) )
!s = nf90_put_var(nc, idVAR  ,  ciccio )
!s = nf90_put_var(nc, idVAR  ,  real(Ens_copy_in,4) )
        if (EnsAveDouble) then
            s = nf90_put_var(nc, idVAR  ,  Ens_copy_in )
        else
            s = nf90_put_var(nc, idVAR  ,  Ens_copy_in_real )
        end if

       call handle_err1(s,counter,fileNetCDF)
!write(*,*) 'prima di chiudere. t: ', mpi_wtime()-tempo

        s =nf90_close(nc)
       call handle_err1(s,counter,fileNetCDF)
!write(*,*) 'chiuso. t: ', mpi_wtime()-tempo
!deallocate(ciccio)

    END SUBROUTINE
       
    SUBROUTINE Ens_WRITE_AVE_BKP(fileNetCDF, VAR,datefrom, dateTo,M, elapsed_time, deflate, deflate_level)
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
       !double precision, allocatable,dimension(:,:,:) :: copy_in
       integer istart,iend

       integer s, nc, counter
       integer timid, depid, yid, xid
       integer idvartime,idgdept,idphit,idlamt,idVAR
       integer shuffle
       integer indexi
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

       !call switch_index_double(M,copy_in,jpiglo,jpjglo,jpk)
        do indexi=1, jpjglo
            Ens_copy_in(:,indexi,:)=transpose(M(:,indexi,:))
        end do
        
       s = nf90_put_var(nc, idVAR  , Ens_copy_in  )

       call handle_err1(s,counter,fileNetCDF)


        s =nf90_close(nc)

    END SUBROUTINE

      SUBROUTINE Ens_trcrst(restart_prefix, ave_freq_1_prefix, ave_freq_2_prefix)
!---------------------------------------------------------------------
!
!                       ROUTINE trcrst
!                     ******************
!
!  PURPOSE :
!  ---------
!     READ files for restart for passive tracer
!
!----------------------------------------------------------------------


       USE calendar
       USE myalloc
       USE TIME_MANAGER
       USE IO_MEM , ONLY : elapsed_time_1, elapsed_time_2, existFilebkp

       IMPLICIT NONE
       
        CHARACTER(len=*), intent(in) :: restart_prefix
        CHARACTER(len=*), intent(in), target :: ave_freq_1_prefix, ave_freq_2_prefix

!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER jn, jn_high
      CHARACTER(LEN=100) filename
      CHARACTER(LEN=100) bkpname
      logical existFile
      logical bkp1hasbeenread,bkp2hasbeenread
      CHARACTER(LEN=:), pointer :: dirct

! 0. initialisations
       bkpname  = 'AVE_FREQ_1/ave.20111231-15:30:00.N1p.nc.bkp'
       filename = 'RESTARTS/RST.20111231-15:30:00.N1p.nc'
       bkp1hasbeenread = .false.
       bkp2hasbeenread = .false.


      IF(lwp) THEN
          WRITE(numout,*) ' '
          WRITE(numout,*) ' *** trcrst beginning of restart for'
          WRITE(numout,*) ' passive tracer'
          WRITE(numout,*) '   with the time TimeStepStart : ',TimeStepStart
          WRITE(numout,*) ' '
      ENDIF

      jn_high=0
      DO jn=1, jptra  ! global loop on tracers to read restart


         filename = restart_prefix//'.'//DateStart//'.'//trim(ctrcnm(jn))// & 
                '.nc'
         CALL readnc_slice_double(trim(filename), 'TRN'//trim(ctrcnm(jn)), trn(:,:,:,jn) )



! ********************   we put initial undef to 0
          trn(:,:,:,jn) = trn(:,:,:,jn) * tmask
      
          trb = trn

!         SECTION AVE backup
!         rsttrn var is used instead of to define another one

          bkpname= ave_freq_2_prefix//'.'//DateStart//'.'//trim(ctrcnm(jn))//'.nc.bkp'

          INQUIRE(FILE=bkpname, EXIST=existFilebkp)

          if (existFilebkp) then
             if (lwp) write(*,*) 'reading ', bkpname
             CALL readnc_slice_double(bkpname,ctrcnm(jn), traIO(:,:,:,jn) )
             if (.not.bkp2hasbeenread) then
               call readnc_scalar_double(bkpname,'elapsed_time',elapsed_time_2)
               call get_att_char(bkpname,'DateStart'  , BKPdatefrom_2)
               bkp2hasbeenread=.true.
             endif
          else
             traIO(:,:,:,jn) = 0.0
          endif



         IF (ctr_hf(jn).eq.1)  THEN
           jn_high = jn_high + 1
           bkpname= ave_freq_1_prefix//'.'//DateStart//'.'//trim(ctrcnm(jn))//'.nc.bkp'
           INQUIRE(FILE=bkpname, EXIST=existFilebkp)

           if (existFilebkp) then
             if (lwp) write(*,*) 'reading ', bkpname
             CALL readnc_slice_double(bkpname,ctrcnm(jn), traIO_HIGH(:,:,:,jn_high) )
                if (.not.bkp1hasbeenread) then
                  call readnc_scalar_double(bkpname,'elapsed_time',elapsed_time_1)
                  call get_att_char(bkpname,'DateStart'  , BKPdatefrom_1)
                  bkp1hasbeenread=.true.
                endif
           else
              traIO_HIGH(:,:,:,jn_high) = 0.0
           endif
         ENDIF


      END DO


! ******************** 3D DIAGNOSTICS  ***********************************
      jn_high=0
      DO jn=1, jptra_dia

          bkpname    = ave_freq_2_prefix//'.'//DateStart//'.'//trim(dianm(jn))//'.nc.bkp'
          INQUIRE(FILE=bkpname, EXIST=existFile)

          if (existFile) then
             if (lwp) write(*,*) 'reading ', bkpname
             CALL readnc_slice_double(bkpname,dianm(jn), tra_DIA_IO(jn,:,:,:) )
             if (.not.bkp2hasbeenread) then
                call readnc_scalar_double(bkpname,'elapsed_time',elapsed_time_2)
                call get_att_char(bkpname,'DateStart'  , BKPdatefrom_2)
                bkp2hasbeenread=.true.
             endif
          else
             tra_DIA_IO(jn,:,:,:) = 0.0
          endif

          IF ((diahf(jn).eq.1).and.(diaWR(jn).eq.1))  THEN
           jn_high = jn_high + 1
           bkpname= ave_freq_1_prefix//'.'//DateStart//'.'//trim(dianm(jn))//'.nc.bkp'
           INQUIRE(FILE=bkpname, EXIST=existFile)
           if (existFile) then
             if (lwp) write(*,*) 'reading ', bkpname
             CALL readnc_slice_double(bkpname,trim(dianm(jn)), tra_DIA_IO_HIGH(jn_high,:,:,:) )
                    if (.not.bkp1hasbeenread) then
                      call readnc_scalar_double(bkpname,'elapsed_time',elapsed_time_1)
                      call get_att_char(bkpname,'DateStart'  , BKPdatefrom_1)
                      bkp1hasbeenread=.true.
                    endif
           else
              tra_DIA_IO_HIGH(jn_high,:,:,:) = 0.0
           endif

          ENDIF
      END DO


! ******************** 2D DIAGNOSTICS  ***********************************

      jn_high=0
      DO jn=1, jptra_dia_2d

          bkpname    = ave_freq_2_prefix//'.'//DateStart//'.'//trim(dianm_2d(jn))//'.nc.bkp'
          INQUIRE(FILE=bkpname, EXIST=existFile)

          if (existFile) then
             if (lwp) write(*,*) 'reading ', bkpname
             CALL readnc_slice_double_2d(bkpname,dianm_2d(jn), tra_DIA_2d_IO(jn,:,:) )
             if (.not.bkp2hasbeenread) then
                call readnc_scalar_double(bkpname,'elapsed_time',elapsed_time_2)
                call get_att_char(bkpname,'DateStart'  , BKPdatefrom_2)
                bkp2hasbeenread=.true.
             endif
          else
             tra_DIA_2d_IO(jn,:,:) = 0.0
          endif

          IF ((diahf_2d(jn).eq.1).and.(diaWR_2d(jn).eq.1))  THEN
           jn_high = jn_high + 1
           bkpname= ave_freq_1_prefix//'.'//DateStart//'.'//trim(dianm_2d(jn))//'.nc.bkp'
           INQUIRE(FILE=bkpname, EXIST=existFile)
           if (existFile) then
             if (lwp) write(*,*) 'reading ', bkpname
             CALL readnc_slice_double_2d(bkpname,trim(dianm_2d(jn)), tra_DIA_2d_IO_HIGH(jn_high,:,:) )
                    if (.not.bkp1hasbeenread) then
                      call readnc_scalar_double(bkpname,'elapsed_time',elapsed_time_1)
                      call get_att_char(bkpname,'DateStart'  , BKPdatefrom_1)
                      bkp1hasbeenread=.true.
                    endif
           else
              tra_DIA_2d_IO_HIGH(jn_high,:,:) = 0.0
           endif

          ENDIF



      END DO




      if(freq_ave_phys .eq. 1) then
        dirct=>ave_freq_1_prefix
      else
        dirct=>ave_freq_2_prefix
      end if
      
      bkpname = dirct//'.'//DateStart//'.'//'vosaline'//'.nc.bkp'
        INQUIRE(FILE=bkpname, EXIST=existFilebkp)
        if (existFilebkp) then
          if(lwp)  write(*,*) 'reading ', bkpname
          call readnc_slice_double(   bkpname,'vosaline',  snIO) 
        else 
          snIO = 0.0
        end if



      bkpname = dirct//'.'//DateStart//'.'//'votemper'//'.nc.bkp'
        INQUIRE(FILE=bkpname, EXIST=existFilebkp)
        if (existFilebkp) then
          if(lwp)  write(*,*) 'reading ', bkpname
          call readnc_slice_double(   bkpname,'votemper',  tnIO)      
        else
          tnIO = 0.0
        end if
     

      bkpname = dirct//'.'//DateStart//'.'//'vozocrtx'//'.nc.bkp'
        INQUIRE(FILE=bkpname, EXIST=existFilebkp)
        if (existFilebkp) then
          if(lwp)  write(*,*) 'reading ', bkpname
          call readnc_slice_double(   bkpname,'vozocrtx',  unIO)  
        else
          unIO = 0.0
        end if


      bkpname = dirct//'.'//DateStart//'.'//'vomecrty'//'.nc.bkp'
        INQUIRE(FILE=bkpname, EXIST=existFilebkp)
        if (existFilebkp) then
          if(lwp)  write(*,*) 'reading ', bkpname
          call readnc_slice_double(   bkpname,'vomecrty',  vnIO)  
        else
          vnIO = 0.0
        end if


      bkpname = dirct//'.'//DateStart//'.'//'vovecrtz'//'.nc.bkp'
        INQUIRE(FILE=bkpname, EXIST=existFilebkp)
        if (existFilebkp) then
          if(lwp)  write(*,*) 'reading ', bkpname
          call readnc_slice_double(   bkpname,'vovecrtz',  wnIO)  
        else
          wnIO = 0.0
        end if


      bkpname = dirct//'.'//DateStart//'.'//'votkeavt'//'.nc.bkp'
        INQUIRE(FILE=bkpname, EXIST=existFilebkp)
        if (existFilebkp) then
          if(lwp)  write(*,*) 'reading ', bkpname
          call readnc_slice_double(   bkpname,'votkeavt',  avtIO)  
        else
          avtIO = 0.0
        end if


      bkpname = dirct//'.'//DateStart//'.'//'e3t'//'.nc.bkp'
        INQUIRE(FILE=bkpname, EXIST=existFilebkp)
        if (existFilebkp) then
          if(lwp)  write(*,*) 'reading ', bkpname
          call readnc_slice_double(   bkpname,'e3t',  e3tIO)  
        else
          e3tIO = 0.0
        end if


      bkpname = dirct//'.'//DateStart//'.'//'soshfldo'//'.nc.bkp'
        INQUIRE(FILE=bkpname, EXIST=existFilebkp)
        if (existFilebkp) then
          if(lwp)  write(*,*) 'reading ', bkpname
          call readnc_slice_double_2d(   bkpname,'soshfldo', qsrIO)  
        else
          qsrIO = 0.0
        end if


      bkpname = dirct//'.'//DateStart//'.'//'sowindsp'//'.nc.bkp'
        INQUIRE(FILE=bkpname, EXIST=existFilebkp)
        if (existFilebkp) then
          if(lwp)  write(*,*) 'reading ', bkpname
          call readnc_slice_double_2d(   bkpname,'sowindsp',  vatmIO)  
        else
          vatmIO = 0.0
        end if


      bkpname = dirct//'.'//DateStart//'.'//'sowaflcd'//'.nc.bkp'
        INQUIRE(FILE=bkpname, EXIST=existFilebkp)
        if (existFilebkp) then
          if(lwp)  write(*,*) 'reading ', bkpname
          call readnc_slice_double_2d(   bkpname,'sowaflcd',  empIO)  
        else
          empIO = 0.0
        end if




      if (lwp) write(*,*) 'trcrst elapsed_time_1 = ', elapsed_time_1, &
                ' elapsed_time_2 = ', elapsed_time_2
      IsStartBackup_1 = bkp1hasbeenread
      IsStartBackup_2 = bkp2hasbeenread


      END SUBROUTINE
      
    SUBROUTINE Ens_trcdia(datemean, datefrom, dateend, FREQ_GROUP)
        use mpi
        
        USE myalloc, &
            only: lwp, trcdiaparttime, trcdiatottime
        use ogstm_mpi_module, &
            only: glcomm
        use Ens_Mem, &
            only: EnsRankZero, &
                EnsDebug, EnsRank, EnsSize, &
                EnsSaveEachAve, EnsSaveMeanAve, &
                Ens_ave_freq_1_prefix, Ens_ave_freq_1_ens_prefix, Ens_ave_freq_2_prefix, Ens_ave_freq_2_ens_prefix, &
                Ens_flux_prefix, Ens_flux_ens_prefix
        use Ens_Utilities, &
            only: int2str, Ens_ReduceMean
        use Ens_Custom, &
            only: Ens_Ave2DA, Ens_DA2Ave

        CHARACTER(LEN=17), INTENT(IN) :: datemean, datefrom, dateend
        INTEGER, INTENT(IN) :: FREQ_GROUP
        
        !local variables
        logical IsBackup
        integer ierr

        if (EnsDebug>0) then
            call mpi_barrier(glcomm, ierr)
            if (lwp) write(*,*) 'Starting Ens_trcdia, group ', FREQ_GROUP
        end if
        trcdiaparttime = MPI_WTIME() ! cronometer-start
        
        IsBackup =  (datemean.eq.dateend)
        
        if(freq_ave_phys.eq.FREQ_GROUP ) then

            if (freq_ave_phys==1)then
                tra_PHYS_2d_IO_high(:,:,1) = vatmIO
                tra_PHYS_2d_IO_high(:,:,2) = empIO
                tra_PHYS_2d_IO_high(:,:,3) = qsrIO
                
                tra_PHYS_IO_high(:,:,:,1) = snIO
                tra_PHYS_IO_high(:,:,:,2) = tnIO
                tra_PHYS_IO_high(:,:,:,3) = wnIO
                tra_PHYS_IO_high(:,:,:,4) = avtIO
                tra_PHYS_IO_high(:,:,:,5) = e3tIO
                tra_PHYS_IO_high(:,:,:,6) = unIO
                tra_PHYS_IO_high(:,:,:,7) = vnIO
            END IF
            
            IF (freq_ave_phys==2)then
                tra_PHYS_2d_IO(:,:,1) = vatmIO
                tra_PHYS_2d_IO(:,:,2) = empIO
                tra_PHYS_2d_IO(:,:,3) = qsrIO
                
                tra_PHYS_IO(:,:,:,1) = snIO
                tra_PHYS_IO(:,:,:,2) = tnIO
                tra_PHYS_IO(:,:,:,3) = wnIO
                tra_PHYS_IO(:,:,:,4) = avtIO
                tra_PHYS_IO(:,:,:,5) = e3tIO
                tra_PHYS_IO(:,:,:,6) = unIO
                tra_PHYS_IO(:,:,:,7) = vnIO
            end if
            
        end if
        
        
        !   writes ave files for tracer concentration
        
        if (EnsSize>1.and.(IsBackup.or.EnsSaveEachAve)) then
            if (EnsDebug>0) then
                call mpi_barrier(glcomm, ierr)
                if (lwp) write(*,*) 'Saving each ave. Time: ', MPI_WTIME()-trcdiaparttime
            end if
            CALL Ens_trcdit(trim(Ens_ave_freq_1_ens_prefix)//int2str(EnsRank,3), & 
                trim(Ens_ave_freq_2_ens_prefix)//int2str(EnsRank,3), &
                datemean, datefrom, dateend,FREQ_GROUP)
            if (EnsDebug>0) then
                call mpi_barrier(glcomm, ierr)
                if (lwp) write(*,*) 'Ens_trcdit done. Time: ', MPI_WTIME()-trcdiaparttime
            end if
            CALL Ens_diadump(trim(Ens_ave_freq_1_ens_prefix)//int2str(EnsRank,3), &
                trim(Ens_ave_freq_2_ens_prefix)//int2str(EnsRank,3), &
                datemean, datefrom, dateend,FREQ_GROUP)
            if (EnsDebug>0) then
                call mpi_barrier(glcomm, ierr)
                if (lwp) write(*,*) 'Ens_diadump done. Time: ', MPI_WTIME()-trcdiaparttime
            end if
            CALL Ens_fluxdump(trim(Ens_flux_ens_prefix)//int2str(EnsRank,3), &
                datemean, datefrom, dateend,FREQ_GROUP)
            if (EnsDebug>0) then
                call mpi_barrier(glcomm, ierr)
                if (lwp) write(*,*) 'Ens_fluxdump done. Time: ', MPI_WTIME()-trcdiaparttime
            end if
        end if
        
        if (EnsSize==1.or.((.not.IsBackup).and.EnsSaveMeanAve)) then
            
            if (EnsSize>1)then
                if (EnsDebug>0) then
                    call mpi_barrier(glcomm, ierr)
                    if (lwp) write(*,*) 'Saving ave mean. Time: ', MPI_WTIME()-trcdiaparttime
                end if
            
                call Ens_Ave2DA
                if (EnsDebug>0) then
                    call mpi_barrier(glcomm, ierr)
                    if (lwp) write(*,*) 'Ave transformed. Time: ', MPI_WTIME()-trcdiaparttime
                end if
                
                call Ens_ReduceMean(win_traIO, n_traIO, gl_traIO)
                call Ens_ReduceMean(win_traIO_HIGH, n_traIO_HIGH, gl_traIO_HIGH)
                
                call Ens_ReduceMean(win_tra_DIA_IO, n_tra_DIA_IO, gl_tra_DIA_IO)
                call Ens_ReduceMean(win_tra_DIA_IO_HIGH, n_tra_DIA_IO_HIGH, gl_tra_DIA_IO_HIGH)
                call Ens_ReduceMean(win_tra_DIA_2d_IO, n_tra_DIA_2d_IO, gl_tra_DIA_2d_IO)
                call Ens_ReduceMean(win_tra_DIA_2d_IO_HIGH, n_tra_DIA_2d_IO_HIGH, gl_tra_DIA_2d_IO_HIGH)
                
                if(freq_ave_phys.eq.FREQ_GROUP ) then
                    if (freq_ave_phys==1) then
                        call Ens_ReduceMean(win_tra_PHYS_IO_HIGH, n_tra_PHYS_IO_HIGH, gl_tra_PHYS_IO_HIGH)
                        call Ens_ReduceMean(win_tra_PHYS_2d_IO_HIGH, n_tra_PHYS_2d_IO_HIGH, gl_tra_PHYS_2d_IO_HIGH)
                    END IF
                    IF (freq_ave_phys==2)then
                        call Ens_ReduceMean(win_tra_PHYS_IO, n_tra_PHYS_IO, gl_tra_PHYS_IO)
                        call Ens_ReduceMean(win_tra_PHYS_2d_IO, n_tra_PHYS_2d_IO, gl_tra_PHYS_2d_IO)
                    END IF
                end if
                
                IF (Fsize.NE.0 ) call Ens_ReduceMean(win_diaflx, n_diaflx, gl_diaflx)
                if (EnsDebug>0) then
                    call mpi_barrier(glcomm, ierr)
                    if (lwp) write(*,*) 'Ave reduced. Time: ', MPI_WTIME()-trcdiaparttime
                end if
                
                call Ens_DA2Ave
                if (EnsDebug>0) then
                    call mpi_barrier(glcomm, ierr)
                    if (lwp) write(*,*) 'Transformed back to ave. Time: ', MPI_WTIME()-trcdiaparttime
                end if
            end if
        
            if (EnsRank==EnsRankZero) CALL Ens_trcdit(trim(Ens_ave_freq_1_prefix), trim(Ens_ave_freq_2_prefix), &
                datemean, datefrom, dateend,FREQ_GROUP)
            if (EnsDebug>0) then
                call mpi_barrier(glcomm, ierr)
                if (lwp) write(*,*) 'Ens_trcdit done. Time: ', MPI_WTIME()-trcdiaparttime
            end if
            if (EnsRank==EnsRankZero) CALL Ens_diadump(trim(Ens_ave_freq_1_prefix), trim(Ens_ave_freq_2_prefix), &
                datemean, datefrom, dateend,FREQ_GROUP)
            if (EnsDebug>0) then
                call mpi_barrier(glcomm, ierr)
                if (lwp) write(*,*) 'Ens_diadump done. Time: ', MPI_WTIME()-trcdiaparttime
            end if
            if (EnsRank==EnsRankZero) CALL Ens_fluxdump(trim(Ens_flux_prefix), datemean, datefrom, dateend,FREQ_GROUP)
            if (EnsDebug>0) then
                call mpi_barrier(glcomm, ierr)
                if (lwp) write(*,*) 'Ens_fluxdump done. Time: ', MPI_WTIME()-trcdiaparttime
            end if
        end if
        
        if (.not.IsBackup) then
            call reset_ave(FREQ_GROUP)
        end if


        trcdiaparttime =   MPI_WTIME() - trcdiaparttime  ! cronometer-stop
        trcdiatottime  = trcdiatottime + trcdiaparttime

    END SUBROUTINE 
    
    SUBROUTINE Ens_trcdit(ave_freq_1_prefix, ave_freq_2_prefix, datemean,datefrom,dateTo,FREQ_GROUP)
        !---------------------------------------------------------------------
        !
        !                       ROUTINE trcdit
        !
        !                     ******************
        !  gcoidess develop
        !
        !  Purpose :
        !  ---------
        !     Standard output of passive tracer : concentration fields



        USE calendar
        USE myalloc
        USE IO_mem
        USE FN_mem
        USE TIME_MANAGER
        use mpi
        USE ogstm_mpi_module
        USE MPI_GATHER_INFO

        USE MATRIX_VARS
        USE NODES_MODULE
        USE DTYPE_PROCS_STRING_MODULE

        IMPLICIT NONE
        CHARACTER(LEN=*), INTENT(IN), target :: ave_freq_1_prefix, ave_freq_2_prefix
        CHARACTER(LEN=17), INTENT(IN) :: datemean, dateFrom, dateTo
        INTEGER, INTENT(IN) :: FREQ_GROUP ! 1 = HIGH FREQ, 2 = LOW FREQ #cannot change value inside te value

        ! local declarations
        ! ==================
        
        double precision ::  Miss_val =1.e20

        INTEGER jk,jj,ji,i,j,k
        INTEGER ind, i_contribution, j_contribution
        double precision :: elapsed_time

        !new declarations
        INTEGER counter_var, counter_var_high, new_counter_var, new_counter_var_high, nVARS, jv, ivar, n_dumping_cycles
        INTEGER col_var, row_var, writing_rank
        CHARACTER(len=20), DIMENSION(nodes) :: matrix_row_to_write

        CHARACTER(LEN=100) output_file_nc  ! AVE_FREQ_1/ave.20091231-12:00:00.P1n.nc
        CHARACTER(LEN=20) var
        CHARACTER(LEN=100) bkpname
        CHARACTER(LEN=:), pointer :: DIR
        logical IsBackup

        CHARACTER(LEN=20)  var_to_store

        !INTEGER :: ind_col
        DOUBLE PRECISION :: start_time_trcdit_info,finish_time_trcdit_info, proctime_time_trcdit_info, max_time_trcdit_info
        DOUBLE PRECISION :: gatherv_fin_time,gatherv_init_time,gatherv_delta_time,gatherv_sum_time,gatherv_mean_time
        DOUBLE PRECISION :: writing_rank_fin_time, writing_rank_init_time,writing_rank_delta_time, writing_rank_sum_time


        !----------------------------------------------------------------------
        ! statement functions
        ! ===================


        INTEGER idrank, ierr, istart, jstart, iPe, iPd, jPe, jPd, status(MPI_STATUS_SIZE)
        INTEGER irange, jrange
        INTEGER totistart, totiend, relistart, reliend
        INTEGER totjstart, totjend, reljstart, reljend

        ! ----------------------------------------
        IsBackup =  (datemean.eq.dateTo)
        if (lwp) write(*,*) 'trcdit IsBackup = ',IsBackup, ' group ' ,FREQ_GROUP
        ! ----------------------------------------
        bkpname  = 'AVE_FREQ_1/ave.20111231-15:30:00.N1p.nc.bkp'
        !call mppsync()


        SELECT CASE (FREQ_GROUP)
        CASE (1) 
                elapsed_time=elapsed_time_1
                DIR=>ave_freq_1_prefix
                n_dumping_cycles = matrix_state_1_row
        CASE (2) 
                elapsed_time=elapsed_time_2
                DIR=>ave_freq_2_prefix
                n_dumping_cycles = matrix_state_2_row
        END SELECT

        if (WRITING_RANK_WR)then
                tottrnIO = Miss_val
        endif
        
        !-----------------------
        !starting big loop
        start_time_trcdit_info= MPI_Wtime()

        COUNTER_VAR = 1
        COUNTER_VAR_HIGH = 1

        DUMPING_LOOP: DO jv = 1, n_dumping_cycles

                gatherv_init_time = MPI_Wtime()

                DO ivar = 1 , nodes!number of variables for each round corresponds to the number of nodes

                        writing_rank = writing_procs(ivar)

                        IF (COUNTER_VAR > JPTRA)then
                                EXIT
                        else if (COUNTER_VAR_HIGH > JPTRA_HIGH)then
                                EXIT
                        ELSE

                                if (FREQ_GROUP.eq.2) then
                                        do ji =1 , jpi
                                                i_contribution= jpk*jpj * (ji - 1 )
                                                do jj =1 , jpj
                                                        j_contribution=jpk*(jj-1)
                                                        do jk =1 , jpk
                                                                ind =  jk + j_contribution + i_contribution
                                                                bufftrn   (ind)= traIO( jk,jj,ji,counter_var)
                                                        enddo
                                                enddo
                                        enddo
                                else ! FREQ_GROUP.eq.1
                                        do ji =1 , jpi
                                                i_contribution= jpk*jpj * (ji - 1 )
                                                do jj =1 , jpj
                                                        j_contribution=jpk*(jj-1)
                                                        do jk =1 , jpk
                                                                ind =  jk + j_contribution + i_contribution
                                                                bufftrn   (ind)= traIO_HIGH(jk,jj,ji,counter_var_high)
                                                        enddo
                                                enddo
                                        enddo
                                endif

                                counter_var = counter_var + 1
                                if (FREQ_GROUP.eq.1) counter_var_high = counter_var_high + 1

                                !GATHERV TO THE WRITING RANK
                                CALL MPI_GATHERV(bufftrn, sendcount, MPI_DOUBLE_PRECISION, bufftrn_TOT, jprcv_count, jpdispl_count, MPI_DOUBLE_PRECISION, writing_rank, mycomm, IERR)

                        END IF

                END DO
                gatherv_fin_time = MPI_Wtime()

                gatherv_delta_time = gatherv_fin_time - gatherv_init_time
                CALL MPI_Reduce( gatherv_delta_time, gatherv_sum_time, 1, MPI_DOUBLE, MPI_SUM, 0, mycomm,IERROR)


                ! *************** COLLECTING DATA *****************************

                IF (WRITING_RANK_WR)then

                        writing_rank_init_time = MPI_Wtime()

                        !ind_col = (myrank / n_ranks_per_node)+1

                        if (FREQ_GROUP.eq.2) then
                                var_to_store = matrix_state_2(jv,ind_col)%var_name
                        else
                                var_to_store = matrix_state_1(jv,ind_col)%var_name
                        end if

                        IF (var_to_store == "novars_input")then
                                EXIT
                        ELSE
                                DO idrank = 0,mysize-1

                                        ! ******* WRITING RANK sets indexes of tot matrix where to place buffers of idrank
                                        irange    = iPe_a(idrank+1) - iPd_a(idrank+1) + 1
                                        jrange    = jPe_a(idrank+1) - jPd_a(idrank+1) + 1
                                        totistart = istart_a(idrank+1) + iPd_a(idrank+1) - 1
                                        totiend   = totistart + irange - 1
                                        totjstart = jstart_a(idrank+1) + jPd_a(idrank+1) - 1
                                        totjend   = totjstart + jrange - 1
                                        relistart = 1 + iPd_a(idrank+1) - 1
                                        reliend   = relistart + irange - 1
                                        reljstart = 1 + jPd_a(idrank+1) - 1
                                        reljend   = reljstart + jrange - 1

                                        ! **** ASSEMBLING *** WRITING RANK  puts in tot matrix buffer received by idrank
                                        do ji =totistart,totiend
                                                i_contribution   = jpk*jpj_rec_a(idrank+1)*(ji-1-totistart+ relistart)
                                                        do jj =totjstart,totjend
                                                                j_contribution = jpk*(jj-1-totjstart+ reljstart)
                                                                        do jk =1, jpk
                                                                                ind = jk + j_contribution + i_contribution
                                                                                tottrnIO(jk,jj,ji)= bufftrn_TOT(ind+jpdispl_count(idrank+1))
                                                                        enddo
                                                        enddo
                                        enddo
                                END DO

                                output_file_nc = DIR//'.'//datemean//'.'//trim(var_to_store)//'.nc'
                                bkpname = DIR//'.'//datemean//'.'//trim(var_to_store)//'.nc.bkp'

                                if (IsBackup) then
                                        CALL Ens_WRITE_AVE_BKP(trim(bkpname),var_to_store,datefrom, dateTo,tottrnIO,elapsed_time, deflate_ave, deflate_level_ave)
                                else
                                        CALL Ens_WRITE_AVE(trim(output_file_nc),var_to_store,datefrom, dateTo, tottrnIO, deflate_ave, deflate_level_ave)
                                endif
                        END IF
                        !writing_rank_fin_time = MPI_Wtime()
                        !writing_rank_delta_time = writing_rank_fin_time - writing_rank_init_time
                        !writing_rank_sum_time = writing_rank_delta_time + writing_rank_sum_time
                        !write(*,*)'writingtime', writing_rank_sum_time,'   ',jv, '   ', myrank
                END IF
        END DO DUMPING_LOOP

        finish_time_trcdit_info= MPI_Wtime()
        proctime_time_trcdit_info=finish_time_trcdit_info - start_time_trcdit_info

        CALL MPI_Reduce( proctime_time_trcdit_info,max_time_trcdit_info, 1, MPI_DOUBLE, MPI_MAX, 0, mycomm,IERROR)

        !if(myrank == 0) then
        !        write(*,*) 'TRCDIT TIME is', max_time_trcdit_info
        !end if


    END SUBROUTINE 
    
    SUBROUTINE Ens_diadump(ave_freq_1_prefix, ave_freq_2_prefix, datemean,datefrom,dateTo,FREQ_GROUP)
!     ******************
      USE calendar
      USE myalloc
      USE IO_mem
      USE FN_mem
      USE TIME_MANAGER
      use mpi
      USE ogstm_mpi_module
      USE MPI_GATHER_INFO
      USE MATRIX_VARS
      USE NODES_MODULE
      USE DTYPE_PROCS_STRING_MODULE

      IMPLICIT NONE

        CHARACTER(LEN=*), INTENT(IN), target :: ave_freq_1_prefix, ave_freq_2_prefix
      CHARACTER(LEN=17), INTENT(IN) :: datemean, dateFrom, dateTo
      INTEGER, INTENT(IN) :: FREQ_GROUP
      INTEGER jk,jj,ji
      INTEGER ind, i_contribution, j_contribution
      CHARACTER(LEN=100) bkpname
      CHARACTER(LEN=:), pointer :: DIR
      logical IsBackup
      double precision :: elapsed_time


      CHARACTER(LEN=100) dia_file_nc
      CHARACTER(LEN=100) phys_file_nc
      CHARACTER(LEN=20)  var

      INTEGER idrank, ierr, istart, jstart, iPe, iPd, jPe, jPd, status(MPI_STATUS_SIZE)
      INTEGER irange, jrange
      INTEGER totistart, totiend, relistart, reliend
      INTEGER totjstart, totjend, reljstart, reljend
      double precision ::  Miss_val =1.e20
      INTEGER :: nVars, counter_var_2d, counter_var_high_2d,counter_var_diag, counter_var_diag_high
      INTEGER :: counter_var_phys_2d,counter_var_phys_high_2d,counter_var_phys, counter_var_phys_high
      CHARACTER(LEN=20) ::  var_to_store_diag_2d, var_to_store_diag
      CHARACTER(LEN=20) ::  var_to_store_phys_2d, var_to_store_phys
      INTEGER :: n_dumping_cycles, jv, ivar, writing_rank!, ind_col
      INTEGER :: var_to_send_2D, var_high_to_send_2D
      INTEGER :: var_to_send, var_high_to_send
     ! call mppsync()
! ----------------------------------------
      IsBackup =  (datemean.eq.dateTo)
      if (lwp) write(*,*) 'diadump IsBackup = ',IsBackup, ' group ' ,FREQ_GROUP
! ----------------------------------------

      if (WRITING_RANK_WR) tottrnIO2d = Miss_val
! ******************  DIAGNOSTIC OUTPUT   2D *******************

        
        IF (FREQ_GROUP==1)then
                elapsed_time=elapsed_time_1
                DIR=>ave_freq_1_prefix
                n_dumping_cycles=matrix_diag_2d_1_row
        end if

        if (FREQ_GROUP==2)then
                elapsed_time=elapsed_time_2
                DIR=>ave_freq_2_prefix
                n_dumping_cycles=matrix_diag_2d_2_row
        END IF


        COUNTER_VAR_2d = 1
        COUNTER_VAR_HIGH_2d = 1


        DUMPING_LOOP_2d: DO jv = 1, n_dumping_cycles

                DO ivar = 1 , nodes

                        writing_rank = writing_procs(ivar)

                        
                        IF (COUNTER_VAR_2d > JPTRA_dia_2d_wri)then
                                EXIT
                        else if (COUNTER_VAR_HIGH_2d > JPTRA_dia_2d_HIGH_wri)then
                                EXIT
                        ELSE

                                var_to_send_2D = lowfreq_table_dia_2d_wri(counter_var_2d)
                                !if (FREQ_GROUP.eq.1) var_high_to_send_2D = highfreq_table_dia_2d_wri(counter_var_high_2d)                        

        
                                if (FREQ_GROUP.eq.2) then
                                        do ji =1 , jpi
                                                i_contribution = jpj * (ji-1)
                                                do jj =1 , jpj
                                                        ind = jj + i_contribution
                                                        buffDIA2d (ind)= tra_DIA_2d_IO(var_to_send_2D,jj,ji)
                                                enddo
                                        enddo
                                else
                                        do ji =1 , jpi
                                                i_contribution = jpj * (ji-1)
                                                do jj = 1 , jpj
                                                        ind = jj + i_contribution
                                                        buffDIA2d (ind)= tra_DIA_2d_IO_high(COUNTER_VAR_HIGH_2d,jj,ji)
                                                enddo
                                        enddo
                                endif
                                counter_var_2d = counter_var_2d + 1
                                if (FREQ_GROUP.eq.1) counter_var_high_2d = counter_var_high_2d + 1

                                CALL MPI_GATHERV(buffDIA2d, sendcount_2d,MPI_DOUBLE_PRECISION, buffDIA2d_TOT,jprcv_count_2d,jpdispl_count_2d, MPI_DOUBLE_PRECISION,writing_rank, mycomm, IERR)

                        END IF
                END DO

        !---------------------------------------------------------------------------------------
        !if writitng rank assembling and dumping

                IF (WRITING_RANK_WR)then

                        !ind_col = (myrank / n_ranks_per_node) +1

                        if (FREQ_GROUP.eq.2) then
                                var_to_store_diag_2d = matrix_diag_2d_2(jv,ind_col)%var_name
                        else
                                var_to_store_diag_2d = matrix_diag_2d_1(jv,ind_col)%var_name
                        end if

                        IF (var_to_store_diag_2d == "novars_input")then
                                EXIT
                        ELSE

                                do idrank = 0,mysize-1
                                         irange    = iPe_a(idrank+1) - iPd_a(idrank+1) + 1
                                         jrange    = jPe_a(idrank+1) - jPd_a(idrank+1) + 1
                                         totistart = istart_a(idrank+1) + iPd_a(idrank+1) - 1
                                         totiend   = totistart + irange - 1
                                         totjstart = jstart_a(idrank+1) + jPd_a(idrank+1) - 1
                                         totjend   = totjstart + jrange - 1
                                         relistart = 1 + iPd_a(idrank+1) - 1
                                         reliend   = relistart + irange - 1
                                         reljstart = 1 + jPd_a(idrank+1) - 1
                                         reljend   = reljstart + jrange - 1
                                         do ji =totistart,totiend ! only 2d vars
                                                i_contribution = jpj_rec_a(idrank+1)*(ji-totistart+relistart -1)
                                                do jj =totjstart,totjend
                                                        ind = jj-totjstart+ reljstart +i_contribution
                                                        tottrnIO2d (jj,ji)=buffDIA2d_TOT(ind+jpdispl_count_2d(idrank+1))
                                                enddo
                                         enddo
                                enddo
                                !if (FREQ_GROUP.eq.2)write(*,*) 'CHECK ', var_to_store_diag_2d,var_to_send_2d
                                !if (FREQ_GROUP.eq.1) write(*,*)'CHECK_h', var_to_store_diag_2d, COUNTER_VAR_HIGH_2d
                                bkpname     = DIR//'.'//datemean//'.'//trim(var_to_store_diag_2d)//'.nc.bkp'
                                dia_file_nc = DIR//'.'//datemean//'.'//trim(var_to_store_diag_2d)//'.nc'

                                if (IsBackup) then
                                        CALL WRITE_AVE_2d_BKP(trim(bkpname),var_to_store_diag_2d,datefrom, dateTo,tottrnIO2d, elapsed_time)

                                else
                                        d2f2d = REAL(tottrnIO2d(:,:),4)
                                        CALL WRITE_AVE_2d(trim(dia_file_nc),var_to_store_diag_2d,datefrom,dateTo, d2f2d)

                                endif
                        end if
                END IF
        END DO DUMPING_LOOP_2d


!------------------------------------------------------------------------------------------------

        if (WRITING_RANK_WR) tottrnIO = Miss_val
! ! ******************  3D DIAGNOSTIC OUTPUT   *******************
        
        IF (FREQ_GROUP==1)then
                !elapsed_time=elapsed_time_1
                !DIR='AVE_FREQ_1/'
                !matrix_col=nodes
                n_dumping_cycles=matrix_diag_1_row
        end if

        if (FREQ_GROUP==2)then
                !elapsed_time=elapsed_time_2
                !DIR='AVE_FREQ_2/'
                !matrix_col=nodes
                n_dumping_cycles=matrix_diag_2_row
        END IF


        COUNTER_VAR_diag = 1
        COUNTER_VAR_diag_HIGH = 1


        DUMPING_LOOP_3d: DO jv = 1, n_dumping_cycles

                DO ivar = 1 , nodes

                        writing_rank = writing_procs(ivar)


                        IF (COUNTER_VAR_diag > JPTRA_dia_wri)then
                                EXIT
                        else if (COUNTER_VAR_diag_HIGH > JPTRA_dia_HIGH_wri)then
                                EXIT
                        ELSE
                                var_to_send = lowfreq_table_dia_wri(counter_var_diag)
                                !if (FREQ_GROUP.eq.1) var_high_to_send = highfreq_table_dia_wri(counter_var_diag_high)

                                if (FREQ_GROUP.eq.2) then
                                        do ji =1, jpi
                                                i_contribution= jpk*jpj * (ji - 1 )
                                                do jj =1 , jpj
                                                        j_contribution=jpk*(jj-1)
                                                        do jk =1 , jpk
                                                                ind = jk + j_contribution + i_contribution
                                                                buffDIA(ind) = tra_DIA_IO(var_to_send, jk,jj,ji)
                                                        enddo
                                                enddo
                                        enddo
                                else
                                        do ji =1, jpi
                                                i_contribution= jpk*jpj * (ji - 1 )
                                                do jj =1 , jpj
                                                        j_contribution=jpk*(jj-1)
                                                        do jk =1 , jpk
                                                                ind = jk + j_contribution + i_contribution
                                                                buffDIA(ind) = tra_DIA_IO_HIGH(COUNTER_VAR_diag_HIGH, jk,jj,ji)
                                                        enddo
                                                enddo
                                        enddo
                                end if
                                !if (FREQ_GROUP.eq.1)write(*,*)'CHECK_h_before', COUNTER_VAR_diag_HIGH
                                counter_var_diag = counter_var_diag + 1
                                if (FREQ_GROUP.eq.1) counter_var_diag_high = counter_var_diag_high + 1

                                !GATHERV TO THE WRITING RANK

                                CALL MPI_GATHERV(buffDIA, sendcount,MPI_DOUBLE_PRECISION, buffDIA_TOT,jprcv_count, jpdispl_count,MPI_DOUBLE_PRECISION, writing_rank,mycomm, IERR)
                        END IF
                END DO

! *********** START WRITING **************************
                IF (WRITING_RANK_WR)then

                        !ind_col = (myrank / n_ranks_per_node)+1

                        if (FREQ_GROUP.eq.2) then
                                var_to_store_diag = matrix_diag_2(jv,ind_col)%var_name
                        else
                                var_to_store_diag = matrix_diag_1(jv,ind_col)%var_name
                        end if

                        IF (var_to_store_diag == "novars_input")then
                                EXIT
                        ELSE

                                do idrank = 0,mysize-1
                                        irange    = iPe_a(idrank+1) - iPd_a(idrank+1) + 1
                                        jrange    = jPe_a(idrank+1) - jPd_a(idrank+1) + 1
                                        totistart = istart_a(idrank+1) + iPd_a(idrank+1) - 1
                                        totiend   = totistart + irange - 1
                                        totjstart = jstart_a(idrank+1) + jPd_a(idrank+1) - 1
                                        totjend   = totjstart + jrange - 1
                                        relistart = 1 + iPd_a(idrank+1) - 1
                                        reliend   = relistart + irange - 1
                                        reljstart = 1 + jPd_a(idrank+1) - 1
                                        reljend   = reljstart + jrange - 1

                                        do ji =totistart,totiend ! 3d vars
                                                i_contribution = jpk*jpj_rec_a(idrank+1)*(ji-1 -totistart+relistart )
                                                do jj =totjstart,totjend
                                                        j_contribution = jpk*(jj-totjstart+ reljstart-1)
                                                        do jk =1 , jpk
                                                                ind = jk + j_contribution + i_contribution
                                                                tottrnIO(jk,jj,ji)= buffDIA_TOT (ind+jpdispl_count(idrank+1))
                                                        enddo
                                                enddo
                                        enddo
                                enddo
                                !if (FREQ_GROUP.eq.2)write(*,*) 'CHECK ', var_to_store_diag, var_to_send
                                !if (FREQ_GROUP.eq.1)write(*,*) 'CHECK_h', var_to_store_diag, COUNTER_VAR_diag_HIGH
                                bkpname     = DIR//'.'//datemean//'.'//trim(var_to_store_diag)//'.nc.bkp'
                                dia_file_nc = DIR//'.'//datemean//'.'//trim(var_to_store_diag)//'.nc'
              
                                if (IsBackup) then
                                        CALL Ens_WRITE_AVE_BKP(trim(bkpname),var_to_store_diag,datefrom, dateTo,tottrnIO,elapsed_time,deflate_ave, deflate_level_ave)
                                else
                                        CALL Ens_WRITE_AVE(trim(dia_file_nc),var_to_store_diag,datefrom,dateTo, tottrnIO,deflate_ave, deflate_level_ave)
                                endif


                        END IF
                END IF
        END DO DUMPING_LOOP_3d

        
!-------------------------------------------------
        ! ****************** PHYSC OUTPUT   2D *******************


        if(freq_ave_phys.eq.FREQ_GROUP ) then

        IF (freq_ave_phys==1)then        
                !elapsed_time=elapsed_time_1
                !DIR='AVE_FREQ_1/'
                n_dumping_cycles=matrix_phys_2d_1_row
        end if

        if (freq_ave_phys==2)then
                !elapsed_time=elapsed_time_2
                !DIR='AVE_FREQ_2/'
                n_dumping_cycles=matrix_phys_2d_2_row
        END IF


        COUNTER_VAR_phys_2d = 1
        COUNTER_VAR_phys_HIGH_2d = 1


        DUMPING_LOOP_2d_phys: DO jv = 1, n_dumping_cycles

                DO ivar = 1 , nodes

                        writing_rank = writing_procs(ivar)

                        !write(*,*)'phys 2d wri number' ,JPTRA_phys_2d_HIGH_wri
                        if (freq_ave_phys==0) then
                            exit
                        else IF (freq_ave_phys==2 .and. COUNTER_VAR_phys_2d > JPTRA_phys_2d_wri)then
                                EXIT
                        else if (freq_ave_phys==1 .and. COUNTER_VAR_phys_HIGH_2d > JPTRA_phys_2d_HIGH_wri)then
                                EXIT
                        ELSE
                                if(freq_ave_phys==2) then
                                        var_to_send_2D = lowfreq_table_phys_2d_wri(counter_var_phys_2d)
                                else
                                        var_to_send_2D = highfreq_table_phys_2d_wri(counter_var_phys_high_2d)
                                        !write(*,*) 'var to send 2d is', var_to_send_2D
                                end if

                                if (freq_ave_phys.eq.2) then
                                        do ji =1 , jpi
                                                i_contribution = jpj * (ji-1)
                                                do jj =1 , jpj
                                                        ind = jj +i_contribution
                                                        buffPHYS2d (ind)=tra_PHYS_2d_IO(jj,ji,var_to_send_2D)
                                                enddo
                                        enddo
                                else
                                        do ji =1 , jpi
                                                i_contribution = jpj * (ji-1)
                                                do jj = 1 , jpj
                                                        ind = jj +i_contribution
                                                        buffPHYS2d (ind)=tra_PHYS_2d_IO_high(jj,ji,var_to_send_2D)
                                                enddo
                                        enddo
                                        !write(*,*) 'valeu in buffer is', tra_PHYS_2d_IO_high(var_to_send_2D,5,5)
                                        !write(*,*) &
                                     !'valeu in buffer is,second print',tra_PHYS_2d_IO_high(3,5,5)
                                endif
                                counter_var_phys_2d = counter_var_phys_2d + 1
                                if (freq_ave_phys.eq.1) counter_var_phys_high_2d = counter_var_phys_high_2d + 1

                                CALL MPI_GATHERV(buffPHYS2d,sendcount_2d,MPI_DOUBLE_PRECISION,buffPHYS2d_TOT,jprcv_count_2d,jpdispl_count_2d,MPI_DOUBLE_PRECISION,writing_rank, mycomm, IERR)

                        END IF
                END DO

        !---------------------------------------------------------------------------------------
        !if writitng rank assembling and dumping

                IF (WRITING_RANK_WR)then

                        !ind_col = (myrank / n_ranks_per_node) +1

                        if (freq_ave_phys.eq.2) then
                                var_to_store_phys_2d = matrix_phys_2d_2(jv,ind_col)%var_name
                        else
                                var_to_store_phys_2d = matrix_phys_2d_1(jv,ind_col)%var_name
                        end if

                        IF (var_to_store_phys_2d == "novars_input")then
                                EXIT
                        ELSE

                                do idrank = 0,mysize-1
                                         irange    = iPe_a(idrank+1) - iPd_a(idrank+1) + 1
                                         jrange    = jPe_a(idrank+1) - jPd_a(idrank+1) + 1
                                         totistart = istart_a(idrank+1) + iPd_a(idrank+1) - 1
                                         totiend   = totistart + irange - 1
                                         totjstart = jstart_a(idrank+1) + jPd_a(idrank+1) - 1
                                         totjend   = totjstart + jrange - 1
                                         relistart = 1 + iPd_a(idrank+1) - 1
                                         reliend   = relistart + irange - 1
                                         reljstart = 1 + jPd_a(idrank+1) - 1
                                         reljend   = reljstart + jrange - 1
                                         do ji =totistart,totiend ! only 2d vars
                                                i_contribution = jpj_rec_a(idrank+1)*(ji-totistart+relistart -1)
                                                do jj =totjstart,totjend
                                                        ind = jj-totjstart+ reljstart +i_contribution
                                                        tottrnIO2d(jj,ji)=buffPHYS2d_TOT(ind+jpdispl_count_2d(idrank+1))
                                                enddo
                                         enddo
                                enddo
                                bkpname     =DIR//'.'//datemean//'.'//trim(var_to_store_phys_2d)//'.nc.bkp'
                                phys_file_nc =DIR//'.'//datemean//'.'//trim(var_to_store_phys_2d)//'.nc'

                                if (IsBackup) then
                                        CALL WRITE_AVE_2d_BKP(trim(bkpname),var_to_store_phys_2d,datefrom, dateTo,tottrnIO2d, elapsed_time)

                                else
                                        d2f2d = REAL(tottrnIO2d(:,:),4)
                                        CALL WRITE_AVE_2d(trim(phys_file_nc),var_to_store_phys_2d,datefrom,dateTo, d2f2d)

                                endif
                        end if
                END IF
        END DO DUMPING_LOOP_2d_phys


        if (WRITING_RANK_WR) tottrnIO = Miss_val
!-------------------------------------------------------------------------


       ! ! ******************  3D PHYS OUTPUT   *******************

        IF (freq_ave_phys==1)then
                !elapsed_time=elapsed_time_1
                !DIR='AVE_FREQ_1/'
                !matrix_col=nodes
                n_dumping_cycles=matrix_phys_1_row
        end if

        if (freq_ave_phys==2)then
                !elapsed_time=elapsed_time_2
                !DIR='AVE_FREQ_2/'
                !matrix_col=nodes
                n_dumping_cycles=matrix_phys_2_row
        END IF


        COUNTER_VAR_phys = 1
        COUNTER_VAR_phys_HIGH = 1


        DUMPING_LOOP_3d_phys: DO jv = 1, n_dumping_cycles

                DO ivar = 1 , nodes

                        writing_rank = writing_procs(ivar)
                        if (freq_ave_phys==0) then
                            exit
                        else IF (freq_ave_phys==2 .and. COUNTER_VAR_phys > JPTRA_phys_wri)then
                                EXIT
                        else if (freq_ave_phys==1 .and. COUNTER_VAR_phys_HIGH > JPTRA_phys_HIGH_wri)then
                                EXIT
                        ELSE
                       
                                if(freq_ave_phys==2) then
                                        var_to_send = lowfreq_table_phys_wri(counter_var_phys)
                                else
                                        var_to_send = highfreq_table_phys_wri(counter_var_phys_high)
                                end if

                                if (freq_ave_phys.eq.2) then
                                        do ji =1, jpi
                                                i_contribution= jpk*jpj * (ji -1 )
                                                do jj =1 , jpj
                                                        j_contribution=jpk*(jj-1)
                                                        do jk =1 , jpk
                                                                ind = jk + j_contribution + i_contribution
                                                                buffPHYS(ind) = tra_PHYS_IO(jk,jj,ji,var_to_send)
                                                        enddo
                                                enddo
                                        enddo
                                else
                                        do ji =1, jpi
                                                i_contribution= jpk*jpj * (ji -1 )
                                                do jj =1 , jpj
                                                        j_contribution=jpk*(jj-1)
                                                        do jk =1 , jpk
                                                                ind = jk + j_contribution + i_contribution
                                                                buffPHYS(ind) = tra_PHYS_IO_HIGH(jk,jj,ji,var_to_send)
                                                        enddo
                                                enddo
                                        enddo
                                end if
                                counter_var_phys = counter_var_phys + 1
                                if (freq_ave_phys.eq.1) counter_var_phys_high =counter_var_phys_high + 1

                                !GATHERV TO THE WRITING RANK

                                CALL MPI_GATHERV(buffPHYS,sendcount,MPI_DOUBLE_PRECISION, buffPHYS_TOT,jprcv_count,jpdispl_count,MPI_DOUBLE_PRECISION, writing_rank,mycomm, IERR)
                        END IF
                END DO

! *********** START WRITING **************************
                IF (WRITING_RANK_WR)then

                        !ind_col = (myrank / n_ranks_per_node)+1

                        if (freq_ave_phys.eq.2) then
                                var_to_store_phys = matrix_phys_2(jv,ind_col)%var_name
                        else
                                var_to_store_phys = matrix_phys_1(jv,ind_col)%var_name
                        end if

                        IF (var_to_store_phys == "novars_input")then
                                EXIT
                        ELSE

                                do idrank = 0,mysize-1
                                        irange    = iPe_a(idrank+1) - iPd_a(idrank+1) + 1
                                        jrange    = jPe_a(idrank+1) - jPd_a(idrank+1) + 1
                                        totistart = istart_a(idrank+1) + iPd_a(idrank+1) - 1
                                        totiend   = totistart + irange - 1
                                        totjstart = jstart_a(idrank+1) + jPd_a(idrank+1) - 1
                                        totjend   = totjstart + jrange - 1
                                        relistart = 1 + iPd_a(idrank+1) - 1
                                        reliend   = relistart + irange - 1
                                        reljstart = 1 + jPd_a(idrank+1) - 1
                                        reljend   = reljstart + jrange - 1

                                        do ji =totistart,totiend ! 3d vars
                                                i_contribution = jpk*jpj_rec_a(idrank+1)*(ji-1 -totistart+relistart )
                                                do jj =totjstart,totjend
                                                        j_contribution = jpk*(jj-totjstart+ reljstart-1)
                                                        do jk =1 , jpk
                                                                ind = jk + j_contribution + i_contribution
                                                                tottrnIO(jk,jj,ji)= buffPHYS_TOT (ind+jpdispl_count(idrank+1))
                                                        enddo
                                                enddo
                                        enddo
                                enddo
                                bkpname     =DIR//'.'//datemean//'.'//trim(var_to_store_phys)//'.nc.bkp'
                                phys_file_nc =DIR//'.'//datemean//'.'//trim(var_to_store_phys)//'.nc'

                                if (IsBackup) then
                                        CALL Ens_WRITE_AVE_BKP(trim(bkpname),var_to_store_phys,datefrom,dateTo,tottrnIO,elapsed_time,deflate_ave, deflate_level_ave)
                                else
                                        CALL Ens_WRITE_AVE(trim(phys_file_nc),var_to_store_phys,datefrom,dateTo, tottrnIO,deflate_ave,deflate_level_ave)
                                endif


                        END IF
                END IF
        END DO DUMPING_LOOP_3d_phys   
        end if

        end SUBROUTINE

      SUBROUTINE Ens_fluxdump(flux_prefix,datemean, datefrom, dateend,FREQ_GROUP)

!     ******************
!     Works only if file Fluxes.nc exists and for FREQ_GROUP=1


      USE DIA_mem
      USE netcdf
      use mpi
      USE IO_mem, only: elapsed_time_1, elapsed_time_2

      IMPLICIT NONE

        CHARACTER(LEN=*), INTENT(IN) :: flux_prefix
      CHARACTER(LEN=17), INTENT(IN) :: datemean, datefrom, dateend
      INTEGER, INTENT(IN) :: FREQ_GROUP ! 1 = HIGH FREQ, 2 = LOW FREQ
      INTEGER jf,js,jn,jn_high, ji,jj,counter,counterV
      INTEGER idrank, ierr, istart,status(MPI_STATUS_SIZE)

      CHARACTER(LEN=100) flux_file


      integer s, nc,INDid
      integer nid, tid
      INTEGER IDS(jptra)
      double precision  Realcounter

      if (.not.existFileFluxes) RETURN
      IF  (freq_flux_dump.eq.1) THEN
           if (FREQ_GROUP.eq.2 )     RETURN
           Realcounter   =    1./elapsed_time_1
      ELSE
           if (FREQ_GROUP.eq.1 )     RETURN
           Realcounter   =    1./elapsed_time_2
      ENDIF

      flx_partTime = MPI_WTIME()

!!!!!!!!!!!!!!!!!!PROCESSOR ZERO COLLECTS ALL FLUXES AND DUMPS  ON NETCDF FILE
      if (myrank ==0) then
!     Rank Zero only creates netcdf file
!     FLUXES/flux.20021216-12:00:00.nc
          flux_file = flux_prefix//'.'//datemean//'.nc'


          ! Just to try without 'or'
          ! s = nf90_create(flux_file, or(nf90_clobber,NF90_HDF5), nc)
          s = nf90_create(trim(flux_file), NF90_HDF5, nc)

          s = nf90_put_att(nc, nf90_global, 'Time_Start'     , datefrom)
          s = nf90_put_att(nc, nf90_global, 'Time___End'     ,  dateend)


          s= nf90_def_dim(nc,'n'           , FsizeGlo,  nid)
          s= nf90_def_dim(nc,'Type'        ,        7,  tid)

          ! ******************   var definition
          s = nf90_def_var(nc,'index' ,   nf90_int, (/nid/),   IndID)

          DO jn=1,jptra
             s = nf90_def_var(nc,ctrcnm(jn) ,   nf90_double, (/tid,nid/),   IDS(jn))
          ENDDO

           s =nf90_enddef(nc)
           counter = 1
           INDflxBuff = 0
           if (Fsize .GT. 0) THEN
               do jf=1,Fsize
                   ji = flx_ridxt(jf,4)
                   jj = flx_ridxt(jf,3)
                   INDflxDUMPZERO(jf)=INDflxDUMP(jf)
                   if ( (ji .EQ. 1) .OR. (ji .EQ. jpi) ) INDflxDUMPZERO(jf) = 0 ! Ghost cell value has index 0
                   if ( (jj .EQ. 1) .OR. (jj .EQ. jpj) ) INDflxDUMPZERO(jf) = 0 ! Ghost cell value has index 0
                   if (INDflxDUMPZERO(jf) .NE. 0) then
                       INDflxDUMPglo(counter) = INDflxDUMPZERO(jf)
                       counter = counter +1
                   endif
               enddo
           endif
           do idrank = 1,mysize-1
               call MPI_RECV(INDflxBuff    , FsizeMax,                 mpi_integer, idrank, 1,mycomm, status, ierr)
               do jf=1,FsizeMax
                   if (INDflxBuff(jf) .NE. 0) then
                       INDflxDUMPglo(counter) = INDflxBuff(jf)
                       counter = counter +1
                   endif
               enddo
           end do
 !     index part
           !************************ now, put var
           s = nf90_put_var(nc, IndID,      INDflxDUMPglo)

           DO jn=1,jptra

               MflxDumpGlo = 0
               counterV=1
               if (Fsize .GT. 0) THEN
                   do jf=1,Fsize
                       if (INDflxDUMPZERO(jf) .NE. 0) then
                           do js=1,7
                               MflxDumpGlo(js,counterV) = diaflx(js,jf,jn)
                           end do
                           counterV = counterV +1
                       endif
                   enddo
               endif
               DO idrank = 1,mysize-1
                   call MPI_RECV(INDflxBuff    , FsizeMax,    mpi_integer, idrank, 2,mycomm, status, ierr)
                   call MPI_RECV(diaflxBuff    , FsizeMax*7,  mpi_real8,   idrank, 3,mycomm, status, ierr)
                   DO jf=1,FsizeMax
                       if (INDflxBuff(jf) .NE. 0) then
                           DO js=1,7
                               MflxDumpGlo(js,counterV) = diaflxBuff(jf,js)
                           END DO
                           counterV = counterV +1
                       endif
                   END DO
               ENDDO ! loop on myrank for each tracers

      !************************ now, put var

          MflxDumpGlo = MflxDumpGlo * Realcounter ! we store average flux
          s = nf90_put_var(nc,IDS(jn) ,  MflxDumpGlo)
!         MflxDumpGlo_FLOAT = real(MflxDumpGlo,4)
!         s = nf90_put_var(nc,IDS(jn) ,  MflxDumpGlo_FLOAT)

          ENDDO ! loop on tracers


          s= nf90_close(nc)

!!!!!!!!!!!!!!!!!!END OF TASK ZERO JOB
      else ! Other ranks than zero
!!!!!!!!!!!!!!!!!!OTHER THAN RANK ZERO JOB STARTS --> SENDING DATA TO PROCESSOR ZERO

           INDflxBuff = 0
           diaflxBuff = 0
           if (Fsize .GT. 0) THEN
               do jf=1,Fsize
                   ji = flx_ridxt(jf,4)
                   jj = flx_ridxt(jf,3)
                   INDflxBuff(jf)=INDflxDUMP(jf)
                   if ( (ji .EQ. 1) .OR. (ji .EQ. jpi) ) INDflxBuff(jf) = 0 ! Ghost cell value has index 0
                   if ( (jj .EQ. 1) .OR. (jj .EQ. jpj) ) INDflxBuff(jf) = 0 ! Ghost cell value has index 0
               enddo
           endif
           call MPI_SEND(INDflxBuff, FsizeMax, mpi_integer, 0, 1, mycomm, ierr)
           DO jn=1,jptra
               IF (Fsize .GT. 0) THEN
                   DO jf=1,Fsize
                       if (INDflxBuff(jf) .NE. 0) then
                           do js=1,7
                               diaflxBuff(jf,js) = diaflx(js,jf,jn)
                           end do
                       endif
                   ENDDO
               ENDIF
               call MPI_SEND(INDflxBuff, FsizeMax, mpi_integer, 0, 2, mycomm, ierr)
               call MPI_SEND(diaflxBuff, FsizeMax*7, mpi_real8, 0, 3, mycomm, ierr)
           END DO ! loop on myrank for each tracers


      endif

!!!!!!!!!!!!!!!!!!END OF COMMUNICATION PART AND FILE DUMP


      flx_partTime = MPI_WTIME() - flx_partTime
      flx_TotTime  = flx_TotTime + flx_partTime

      end SUBROUTINE
      
    Subroutine reset_ave(FREQ_GROUP)
    
        integer, intent(in) :: FREQ_GROUP
                
        if (FREQ_GROUP.eq.2) then
            traIO(:,:,:,:) = 0.0d0
            tra_DIA_2d_IO(:,:,:) = 0.0d0
            tra_DIA_IO(:,:,:,:) = 0.0d0
        else
            traIO_HIGH(:,:,:,:) = 0.0d0
            tra_DIA_2d_IO_HIGH(:,:,:) = 0.0d0
            tra_DIA_IO_HIGH(:,:,:,:) = 0.0d0
        endif

        if ( freq_ave_phys.eq.FREQ_GROUP) then
            snIO     = 0.0d0
            tnIO     = 0.0d0
            vatmIO   = 0.0d0
            empIO    = 0.0d0
            qsrIO    = 0.0d0
            unIO     = 0.0d0
            vnIO     = 0.0d0
            wnIO     = 0.0d0
            avtIO    = 0.0d0
            e3tIO    = 00d0
        endif
        
        IF (Fsize.NE.0 ) diaflx = 0
      
    end SUBROUTINE



end module
