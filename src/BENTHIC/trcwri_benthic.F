      SUBROUTINE trcwri_benthic(datestring)


      USE myalloc
      USE IO_mem
      USE BIOBEN_mem
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

       trcwriparttime = MPI_WTIME() ! cronometer-start

      call mppsync()

      buffben3D     = Miss_val
      bufftrn_ben   = Miss_val
      bufftrb_ben   = Miss_val

      DO jn=1,jptra_b
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


            do jk =1 , jpk_ben
            do jj =1 , jpj
            do ji =1 , jpi
               if (tmask(ji,jj,1).eq.1.0) buffben3D(ji,jj, jk) = trn_ben(ji,jj, jk, jn)
            enddo
            enddo
            enddo
           tottrn_ben (totistart:totiend, totjstart:totjend,:)= buffben3D (relistart:reliend, reljstart:reljend, :)


            do jk =1 , jpk_ben
            do jj =1 , jpj
            do ji =1 , jpi
               if (tmask(ji,jj,1).eq.1.0) buffben3D(ji,jj, jk) = trb_ben(ji,jj, jk, jn)
            enddo
            enddo
            enddo

           tottrb_ben (totistart:totiend, totjstart:totjend,:)= buffben3D  (relistart:reliend, reljstart:reljend, :)




           do idrank = 1, size-1

              call MPI_RECV(jpi_rec , 1,                  mpi_integer, idrank,  1,mpi_comm_world, status, ierr)
              call MPI_RECV(jpj_rec , 1,                  mpi_integer, idrank,  2,mpi_comm_world, status, ierr)
              call MPI_RECV(istart  , 1,                  mpi_integer, idrank,  3,mpi_comm_world, status, ierr)
              call MPI_RECV(jstart  , 1,                  mpi_integer, idrank,  4,mpi_comm_world, status, ierr)
              call MPI_RECV(iPe     , 1,                  mpi_integer, idrank,  5,mpi_comm_world, status, ierr)
              call MPI_RECV(jPe     , 1,                  mpi_integer, idrank,  6,mpi_comm_world, status, ierr)
              call MPI_RECV(iPd     , 1,                  mpi_integer, idrank,  7,mpi_comm_world, status, ierr)
              call MPI_RECV(jPd     , 1,                  mpi_integer, idrank,  8,mpi_comm_world, status, ierr)
              call MPI_RECV(bufftrn_ben,jpi_rec*jpj_rec*jpk_ben, mpi_real8, idrank, 11,mpi_comm_world, status, ierr)
              call MPI_RECV(bufftrb_ben,jpi_rec*jpj_rec*jpk_ben, mpi_real8, idrank, 12,mpi_comm_world, status, ierr)

              write(*,*) "in trcwri_benthic.F: idrank", idrank , "jpi_rec*jpj_rec*jpk_ben",jpi_rec*jpj_rec*jpk_ben 

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


              do jk =1 , jpk_ben
                  do jj =totjstart,totjend
                  ind1=(                   relistart +(jj-totjstart + reljstart-1)*jpi_rec + (jk-1)*jpj_rec*jpi_rec)
                  ind2=((totiend-totistart+relistart)+(jj-totjstart + reljstart-1)*jpi_rec + (jk-1)*jpj_rec*jpi_rec)
                  tottrn_ben(totistart:totiend, jj, jk) =bufftrn_ben(ind1:ind2)
                  tottrb_ben(totistart:totiend, jj, jk) =bufftrb_ben(ind1:ind2)
                  enddo
              enddo


           enddo ! do idrank = 1, size-1



        else !rank != 0


            do jk =1 , jpk_ben
            do jj =1 , jpj
            do ji =1 , jpi
               ind1 = ji + jpi * (jj-1) + jpi * jpj *(jk-1)
               if (tmask(ji,jj,1).eq.1.0) then
                  bufftrn_ben(ind1)= trn_ben(ji,jj, jk, jn)
                  bufftrb_ben(ind1)= trb_ben(ji,jj, jk, jn)
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
            call MPI_SEND(bufftrn_ben  ,jpi*jpj*jpk_ben,  mpi_real8, 0, 11, mpi_comm_world,ierr)
            call MPI_SEND(bufftrb_ben  ,jpi*jpj*jpk_ben,  mpi_real8, 0, 12, mpi_comm_world,ierr)

        endif ! if rank = 0


        if(rank == 0) then

            varname=cbennm(jn)
            filename = 'RESTARTS/RST.'//datestring//'.'//varname//'.nc'

        CALL write_restart_ben(filename,julian)

        endif ! if rank = 0
      END DO ! DO jn=1,jptra_b


       trcwriparttime = MPI_WTIME() - trcwriparttime
       trcwritottime = trcwritottime + trcwriparttime



      END SUBROUTINE trcwri_benthic
