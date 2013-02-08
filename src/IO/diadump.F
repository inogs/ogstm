      SUBROUTINE diadump(datemean,datefrom,dateTo)
!     ******************
      USE calendar
      USE myalloc
      USE IO_mem
      USE FN_mem
      USE TIME_MANAGER
#ifdef key_mpp
      USE myalloc_mpp
#endif

      IMPLICIT NONE


      CHARACTER(LEN=17), INTENT(IN) :: datemean, dateFrom, dateTo

      INTEGER ji, jj, jk, jn
      INTEGER ind

      CHARACTER(LEN=29) forcing_file
      CHARACTER(LEN=32) bkpname
      logical IsBackup


      CHARACTER(LEN=28) dia_file_nc
      CHARACTER(LEN=3)  var

      INTEGER idrank, ierr, istart, jstart, iPe, iPd, jPe, jPd, status(MPI_STATUS_SIZE)
      INTEGER irange, jrange
      INTEGER totistart, totiend, relistart, reliend
      INTEGER totjstart, totjend, reljstart, reljend


      call mppsync()
! ----------------------------------------
      IsBackup =  (datemean.eq.dateTo)
      if (lwp) write(*,*) 'diadump IsBackup = ',IsBackup
! ----------------------------------------
      bkpname  = 'ave.20111231-15:30:00.N1p.nc.bkp'




! ******************  DIAGNOSTIC OUTPUT   *******************
      DO jn =1 , jptra_dia

          if (rank == 0) then                    ! IF LABEL 1


! ******* rank 0 sets indexes of tot matrix where to place its own part

             iPd    = nldi
             iPe    = nlei
             jPd    = nldj
             jPe    = nlej
             istart = nimpp
             jstart = njmpp
             irange    = iPe - iPd + 1
             jrange    = jPe - jPd + 1
             totistart = istart + iPd - 1 ;totiend   = totistart + irange - 1
             totjstart = jstart + jPd - 1 ; totjend   = totjstart + jrange - 1
             relistart = 1 + iPd - 1      ; reliend   = relistart + irange - 1
             reljstart = 1 + jPd - 1      ; reljend   = reljstart + jrange - 1

             tottrnIO (totistart:totiend, totjstart:totjend,:) = tra_ppIO(relistart:reliend, reljstart:reljend, :,jn) ! diagnostic from reaction model

             do idrank = 1, size-1
! **************  rank 0 is receiving from the others their buffer  ****

                call MPI_RECV(jpi_rec    , 1,                 mpi_integer, idrank, 22,mpi_comm_world, status, ierr) !* first info to know where idrank is working
                call MPI_RECV(jpj_rec    , 1,                 mpi_integer, idrank, 23,mpi_comm_world, status, ierr)
                call MPI_RECV(istart     , 1,                 mpi_integer, idrank, 24,mpi_comm_world, status, ierr)
                call MPI_RECV(jstart     , 1,                 mpi_integer, idrank, 25,mpi_comm_world, status, ierr)
                call MPI_RECV(iPe        , 1,                 mpi_integer, idrank, 26,mpi_comm_world, status, ierr)
                call MPI_RECV(jPe        , 1,                 mpi_integer, idrank, 27,mpi_comm_world, status, ierr)
                call MPI_RECV(iPd        , 1,                 mpi_integer, idrank, 28,mpi_comm_world, status, ierr)
                call MPI_RECV(jPd        , 1                 ,mpi_integer, idrank, 29,mpi_comm_world, status, ierr)

                call MPI_RECV(buffpp  ,jpi_rec*jpj_rec*jpk*jptra_dia,mpi_real8,idrank, 30,mpi_comm_world, status, ierr)

! ******* rank 0 sets indexes of tot matrix where to place buffers of idrank
                irange    = iPe - iPd + 1
                jrange    = jPe - jPd + 1
                totistart = istart + iPd - 1 ; totiend   = totistart + irange - 1
                totjstart = jstart + jPd - 1 ; totjend   = totjstart + jrange - 1
                relistart = 1 + iPd - 1      ; reliend   = relistart + irange - 1
                reljstart = 1 + jPd - 1      ; reljend   = reljstart + jrange - 1


                do jk =1 , jpk ! 3d vars
                     do jj =totjstart,totjend
                        do ji =totistart,totiend
                           ind = (ji-totistart+ relistart )+ (jj-totjstart+ reljstart-1)*jpi_rec+(jk-1)*jpj_rec*jpi_rec
                           tottrnIO(ji,jj,jk)= buffpp (ind)
                        enddo
                     enddo
                  enddo


             enddo !idrank = 1, size-1


          else  ! IF LABEL 1,  if(rank == 0)


            do jk =1, jpk
             do jj =1 , jpj
              do ji =1 , jpi
                  ind         =  ji + jpi * (jj-1) + jpi * jpj *(jk-1)
                  buffpp(ind) = tra_ppIO(ji,jj, jk,jn)
              enddo
             enddo
            enddo



              call MPI_SEND(jpi  , 1,mpi_integer, 0, 22, mpi_comm_world,ierr)
              call MPI_SEND(jpj  , 1,mpi_integer, 0, 23, mpi_comm_world,ierr)
              call MPI_SEND(nimpp, 1,mpi_integer, 0, 24, mpi_comm_world,ierr)
              call MPI_SEND(njmpp, 1,mpi_integer, 0, 25, mpi_comm_world,ierr)
              call MPI_SEND(nlei , 1,mpi_integer, 0, 26, mpi_comm_world,ierr)
              call MPI_SEND(nlej , 1,mpi_integer, 0, 27, mpi_comm_world,ierr)
              call MPI_SEND(nldi , 1,mpi_integer, 0, 28, mpi_comm_world,ierr)
              call MPI_SEND(nldj , 1,mpi_integer, 0, 29, mpi_comm_world,ierr)

              call MPI_SEND(buffpp  , jpi*jpj*jpk,mpi_real8, 0, 30, mpi_comm_world,ierr)



      endif ! IF LABEL 1, if(rank == 0)
!************* END COLLECTING DATA  *****************

! *********** START WRITING **************************

      if (rank == 0) then
              var        =  dia(jn)
              bkpname     = 'ave.'//datemean//'.'//var//'.nc.bkp'
              dia_file_nc = 'ave.'//datemean//'.'//var//'.nc'

              if (IsBackup) then
                 CALL WRITE_AVE_BKP(bkpname,datefrom, dateTo,tottrnIO(:,:,:));
                 !CALL IOOGSNC_AVE_BKP(bkpname,var,jpiglo,jpjglo,jpk,lon,lat,dep,datefrom,dateTo,totdiaIO(:,:,:,jn),AC)
              else
                 d2f3d = REAL(tottrnIO(:,:,:),4)
                 CALL WRITE_AVE(dia_file_nc,datefrom,dateTo, d2f3d);
                 !CALL IOOGSNC_AVE_3D2(dia_file_nc, var,jpiglo,jpjglo,jpk,lon,lat,dep,datefrom, dateTo,d2f3d)

              endif


      end if ! IF LABEL 4  if(rank == 0)

      enddo  ! loop in jn

      if (.not.IsBackup) then
          snIO     = 0.
          tnIO     = 0.
          vatmIO   = 0.
          empIO    = 0.
          qsrIO    = 0.
          unIO     = 0.
          bblxIO   = 0.
          vnIO     = 0.
          bblyIO   = 0.
          wnIO     = 0.
          avtIO    = 0.
          tra_ppIO = 0.
      endif




      end SUBROUTINE diadump
