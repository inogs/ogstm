      SUBROUTINE trcwriDA(datestring)

      USE netcdf
      USE myalloc
      USE IO_mem
      USE calendar
      USE TIME_MANAGER
      use mpi
      USE ogstm_mpi_module
      USE DA_mem



      IMPLICIT NONE
      CHARACTER(LEN=17), INTENT(IN) :: datestring



!----------------------------------------------------------------------
! local declarations
! ==================
      double precision ::  Miss_val =1.e20
      INTEGER jk,jj,ji,jn
      INTEGER s, nc, counter
      integer timid, depid, yid, xid, idvar
      double precision julian


      CHARACTER(LEN=45) BeforeName
      CHARACTER(LEN=43) BeforeNameShort

      CHARACTER(LEN=3) varname

      INTEGER idrank, ierr, istart, jstart, iPe, iPd, jPe, jPd, status(MPI_STATUS_SIZE)
      INTEGER irange, jrange
      INTEGER totistart, totiend, relistart, reliend
      INTEGER totjstart, totjend, reljstart, reljend
      INTEGER ind1, i_contribution, j_contribution
      INTEGER SysErr, system



       julian=datestring2sec(datestring)

       if(lwp)write(*,*) 'trcwri DA ------------  myrank =',myrank,' datestring = ',  datestring

       trcwriparttime = MPI_WTIME() ! cronometer-start

      call mppsync()
      if (lwp)  CHLtot = 0.0
      buf     = Miss_val
      bufftrn = Miss_val
      if (lwp) tottrn = Miss_val

       DO jn=1,jptra
        if(.not.isaDAvar(ctrcnm(jn))) CYCLE
        if(myrank == 0) then
           istart = nimpp
           jstart = njmpp
           iPd = nldi
           iPe = nlei
           jPd = nldj
           jPe = nlej
           irange    = iPe - iPd + 1
           jrange    = jPe - jPd + 1
           totistart = istart + iPd - 1
           totiend   = totistart + irange - 1
           totjstart = jstart + jPd - 1
           totjend   = totjstart + jrange - 1
           relistart = 1 + iPd - 1
           reliend   = relistart + irange - 1
           reljstart = 1 + jPd - 1
           reljend   = reljstart + jrange - 1


            do ji =1 , jpi
                  do jj =1 , jpj
                        do jk =1 , jpk
                              if (tmask(jk,jj,ji).eq.1) then
                                    buf(jk,jj,ji) = trn(jk,jj,ji, jn)
                              endif
                        enddo
                  enddo
            enddo
           tottrn  (:,totjstart:totjend,totistart:totiend)= buf(:, reljstart:reljend,relistart:reliend)




      do idrank = 1,mpi_glcomm_size-1
         call MPI_RECV(jpi_rec , 1,                  mpi_integer, idrank,  1,mpi_comm_world, status, ierr)
         call MPI_RECV(jpj_rec , 1,                  mpi_integer, idrank,  2,mpi_comm_world, status, ierr)
         call MPI_RECV(istart  , 1,                  mpi_integer, idrank,  3,mpi_comm_world, status, ierr)
         call MPI_RECV(jstart  , 1,                  mpi_integer, idrank,  4,mpi_comm_world, status, ierr)
         call MPI_RECV(iPe     , 1,                  mpi_integer, idrank,  5,mpi_comm_world, status, ierr)
         call MPI_RECV(jPe     , 1,                  mpi_integer, idrank,  6,mpi_comm_world, status, ierr)
         call MPI_RECV(iPd     , 1,                  mpi_integer, idrank,  7,mpi_comm_world, status, ierr)
         call MPI_RECV(jPd     , 1,                  mpi_integer, idrank,  8,mpi_comm_world, status, ierr)
         call MPI_RECV(bufftrn,   jpi_rec*jpj_rec*jpk, mpi_real8, idrank, 11,mpi_comm_world, status, ierr)



           irange    = iPe - iPd + 1
           jrange    = jPe - jPd + 1
           totistart = istart + iPd - 1
           totiend   = totistart + irange - 1
           totjstart = jstart + jPd - 1
           totjend   = totjstart + jrange - 1
           relistart = 1 + iPd - 1
           reliend   = relistart + irange - 1
           reljstart = 1 + jPd - 1
           reljend   = reljstart + jrange - 1

          do ji =totistart,totiend
            i_contribution   = jpk*jpj_rec*(ji-1-totistart+ relistart)
            do jj =totjstart,totjend
              j_contribution = jpk*(jj-1-totjstart+ reljstart)
              do jk =1, jpk
                  ind1 = jk + j_contribution + i_contribution
                  tottrn(jk,jj,ji)= bufftrn(ind1)
              enddo
            enddo
          enddo

          enddo ! do idrank = 1, size-1
       else !myrank != 0
           do ji =1 , jpi
            i_contribution= jpk*jpj * (ji - 1 )
            do jj =1 , jpj
             j_contribution=jpk*(jj-1)
             do jk =1 , jpk
              ind1 = jk + j_contribution + i_contribution
              if (tmask(jk,jj,ji).eq.1) then
                 bufftrn(ind1)= trn(jk,jj,ji, jn)
              endif
             enddo
            enddo
           enddo
           call MPI_SEND(jpi      , 1         ,mpi_integer, 0,  1, mpi_comm_world,ierr)
           call MPI_SEND(jpj      , 1         ,mpi_integer, 0,  2, mpi_comm_world,ierr)
           call MPI_SEND(nimpp    , 1         ,mpi_integer, 0,  3, mpi_comm_world,ierr)
           call MPI_SEND(njmpp    , 1         ,mpi_integer, 0,  4, mpi_comm_world,ierr)
           call MPI_SEND(nlei     , 1         ,mpi_integer, 0,  5, mpi_comm_world,ierr)
           call MPI_SEND(nlej     , 1         ,mpi_integer, 0,  6, mpi_comm_world,ierr)
           call MPI_SEND(nldi     , 1         ,mpi_integer, 0,  7, mpi_comm_world,ierr)
           call MPI_SEND(nldj     , 1         ,mpi_integer, 0,  8, mpi_comm_world,ierr)
           call MPI_SEND(bufftrn  ,jpk*jpj*jpi,  mpi_real8, 0, 11, mpi_comm_world,ierr)
       endif ! if myrank = 0


        if(myrank == 0) then

            varname=ctrcnm(jn)
            BeforeName = 'DA__FREQ_1/RSTbefore.'//datestring//'.'//varname//'.nc'
            BeforeNameShort = 'DA__FREQ_1/RSTbefore.'//datestring(1:11)//datestring(13:14)//datestring(16:17)//'.'//varname//'.nc'
            do ji=1,jpiglo
              do jj=1,jpjglo
                  do jk=1,jpk
                    tottrnDA(ji,jj,jk) = REAL(tottrn(jk,jj,ji),4)
                  end do
              end do
            end do
            CALL write_BeforeAss(BeforeName, varname)

            !process 0 creates link to the restart file
            !since parallel-netcdf seems to do not 
            !read filenames with colons
            SysErr = system("ln -sf $PWD/"//BeforeName//" "//BeforeNameShort)
            if(SysErr /= 0) call MPI_Abort(MPI_COMM_WORLD, -1, SysErr)
            
            if (isaCHLVAR(varname)) then
              do jk=1,jpk
                do jj=1,jpjglo
                  do ji=1,jpiglo
                    CHLtot(ji,jj,jk) = CHLtot(ji,jj,jk) + tottrnDA(ji,jj,jk)
                  end do
                end do
              end do
            endif
           

        endif ! if myrank = 0
      END DO ! DO jn=1,jptra

      if(myrank == 0) then

        s = nf90_create(CHLSUP_FOR_DA, NF90_CLOBBER, nc)
        s= nf90_def_dim(nc,'lon'           , jpiglo,  xid)
        s= nf90_def_dim(nc,'lat'           , jpjglo,  yid)
        s= nf90_def_dim(nc,'depth'         , 1   ,depid)
        s= nf90_def_dim(nc,'time'  , NF90_UNLIMITED,timid)
        s = nf90_def_var(nc, 'lchlm', nf90_float, (/xid,yid,depid,timid/),  idVAR)

        s = nf90_put_att(nc,idVAR, 'missing_value',4.e+20)

        s =nf90_enddef(nc)

        CHL_SUP = CHLtot(:,:,1)
        s = nf90_put_var(nc, idVAR  ,CHL_SUP ); call handle_err1(s,counter,CHLSUP_FOR_DA)
        s=nf90_close(nc)

      endif

       trcwriparttime = MPI_WTIME() - trcwriparttime
       trcwritottime = trcwritottime + trcwriparttime



      END SUBROUTINE trcwriDA

       !****************************************************************************
       !****************************************************************************
       !****************************************************************************
!      writes tottrnDA as float on NetCDF classic file, because at the moment
!      3d_var uses parallel netcdf that does not work with NetCDF4
       SUBROUTINE write_BeforeAss(fileNetCDF, VAR)

       USE netcdf
       USE myalloc
       USE DA_mem

       IMPLICIT NONE
       CHARACTER*(*), intent(in) :: fileNetCDF
       CHARACTER(LEN=3),intent(in):: VAR

       ! local

       integer s, nc, counter
       integer depid, yid, xid,idN

        s = nf90_create(fileNetCDF, NF90_CLOBBER, nc)

        ! *********** DIMENSIONS ****************
        s= nf90_def_dim(nc,'x'   , jpiglo,  xid)
        s= nf90_def_dim(nc,'y'   , jpjglo,  yid)
        s= nf90_def_dim(nc,'z'   , jpk   ,depid)


        s = nf90_def_var(nc,VAR, nf90_float, (/xid,yid,depid /), idN)
        s = nf90_put_att(nc,idN   , 'missing_value',1.e+20)
        s =nf90_enddef(nc)
        s = nf90_put_var(nc, idN,  tottrnDA); call handle_err1(s,counter,fileNetCDF)
        s =nf90_close(nc)


       END SUBROUTINE write_BeforeAss


