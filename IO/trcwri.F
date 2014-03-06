      SUBROUTINE trcwri(datestring)


      USE myalloc
      USE IO_mem
      USE calendar
      USE TIME_MANAGER
#ifdef key_mpp
      USE myalloc_mpp
#endif

      IMPLICIT NONE
      CHARACTER(LEN=17), INTENT(IN) :: datestring



!----------------------------------------------------------------------
! local declarations
! ==================
      REAL(8) ::  Miss_val =1.e20
      INTEGER ji,jj,jk,jn
      REAL(8) julian


      CHARACTER(LEN=37) filename

      CHARACTER(LEN=3) varname

      INTEGER idrank, ierr, istart, jstart, iPe, iPd, jPe, jPd, status(MPI_STATUS_SIZE)
      INTEGER irange, jrange
      INTEGER totistart, totiend, relistart, reliend
      INTEGER totjstart, totjend, reljstart, reljend
      INTEGER ind1, ind2


       filename = 'RST.20111231-15:30:00.N1p.nc'

       julian=datestring2sec(datestring)

       if(lwp)write(*,*) 'trcwri ------------  rank =',rank,' datestring = ',  datestring

       trcwriparttime = MPI_WTIME() ! F79 cronometer-start

      call mppsync()

      buf     = Miss_val
      bufftrb = Miss_val
      bufftrn = Miss_val

      DO jn=1,jptra
        if(rank == 0) then
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


            do jk =1 , jpk
            do jj =1 , jpj
            do ji =1 , jpi
               if (tmask(ji,jj,jk).eq.1.0) buf(ji,jj, jk) = trn(ji,jj, jk, jn)
            enddo
            enddo
            enddo
           tottrn  (totistart:totiend, totjstart:totjend,:)= buf  (relistart:reliend, reljstart:reljend, :)


            do jk =1 , jpk
            do jj =1 , jpj
            do ji =1 , jpi
               if (tmask(ji,jj,jk).eq.1.0) buf(ji,jj, jk) = trb(ji,jj, jk, jn)
            enddo
            enddo
            enddo

           tottrb  (totistart:totiend, totjstart:totjend,:)= buf  (relistart:reliend, reljstart:reljend, :)




           do idrank = 1, size-1

              call MPI_RECV(jpi_rec , 1,                  mpi_integer, idrank,  1,mpi_comm_world, status, ierr)
              call MPI_RECV(jpj_rec , 1,                  mpi_integer, idrank,  2,mpi_comm_world, status, ierr)
              call MPI_RECV(istart  , 1,                  mpi_integer, idrank,  3,mpi_comm_world, status, ierr)
              call MPI_RECV(jstart  , 1,                  mpi_integer, idrank,  4,mpi_comm_world, status, ierr)
              call MPI_RECV(iPe     , 1,                  mpi_integer, idrank,  5,mpi_comm_world, status, ierr)
              call MPI_RECV(jPe     , 1,                  mpi_integer, idrank,  6,mpi_comm_world, status, ierr)
              call MPI_RECV(iPd     , 1,                  mpi_integer, idrank,  7,mpi_comm_world, status, ierr)
              call MPI_RECV(jPd     , 1,                  mpi_integer, idrank,  8,mpi_comm_world, status, ierr)
              call MPI_RECV(bufftrn,   jpi_rec*jpj_rec*jpk, mpi_real8, idrank, 11,mpi_comm_world, status, ierr)
              call MPI_RECV(bufftrb,   jpi_rec*jpj_rec*jpk, mpi_real8, idrank, 12,mpi_comm_world, status, ierr)


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


              do jk =1 , jpk
                  do jj =totjstart,totjend
                  ind1=(                   relistart +(jj-totjstart + reljstart-1)*jpi_rec + (jk-1)*jpj_rec*jpi_rec)
                  ind2=((totiend-totistart+relistart)+(jj-totjstart + reljstart-1)*jpi_rec + (jk-1)*jpj_rec*jpi_rec)
                  tottrn(totistart:totiend, jj, jk) =bufftrn(ind1:ind2)
                  tottrb(totistart:totiend, jj, jk) =bufftrb(ind1:ind2)
                  enddo
              enddo


           enddo ! do idrank = 1, size-1



        else !rank != 0


            do jk =1 , jpk
            do jj =1 , jpj
            do ji =1 , jpi
               ind1 = ji + jpi * (jj-1) + jpi * jpj *(jk-1)
               if (tmask(ji,jj,jk).eq.1.0) then
                  bufftrn(ind1)= trn(ji,jj, jk, jn)
                  bufftrb(ind1)= trb(ji,jj, jk, jn)
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
            call MPI_SEND(bufftrn  ,jpi*jpj*jpk,  mpi_real8, 0, 11, mpi_comm_world,ierr)
            call MPI_SEND(bufftrb  ,jpi*jpj*jpk,  mpi_real8, 0, 12, mpi_comm_world,ierr)

        endif ! if rank = 0


        if(rank == 0) then

            varname=ctrcnm(jn)
            filename = 'RESTARTS/RST.'//datestring//'.'//varname//'.nc'

        CALL write_restart(filename,julian)

        endif ! if rank = 0
      END DO ! DO jn=1,jptra


!      qui scriveva il date.step

       trcwriparttime = MPI_WTIME() - trcwriparttime
       trcwritottime = trcwritottime + trcwriparttime



      END SUBROUTINE trcwri
