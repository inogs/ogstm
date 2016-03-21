      SUBROUTINE diawri_benthic(datemean,datefrom,dateTo,FREQ_GROUP)
!     ******************
      USE calendar
      USE myalloc
      USE IO_mem
      USE BIOBEN_mem
      USE FN_mem
      USE TIME_MANAGER
#ifdef key_mpp
      USE myalloc_mpp
#endif

      IMPLICIT NONE


      CHARACTER(LEN=17), INTENT(IN) :: datemean, dateFrom, dateTo
      INTEGER, INTENT(IN) :: FREQ_GROUP
      INTEGER ji, jj, jk, jn, jn_high
      INTEGER ind

      CHARACTER(LEN=42) forcing_file
      CHARACTER(LEN=60) bkpname
      CHARACTER(LEN=11) DIR
      logical IsBackup
      integer ave_counter


      CHARACTER(LEN=56) dia_file_nc
      CHARACTER(LEN=20)  var

      INTEGER idrank, ierr, istart, jstart, iPe, iPd, jPe, jPd, status(MPI_STATUS_SIZE)
      INTEGER irange, jrange
      INTEGER totistart, totiend, relistart, reliend
      INTEGER totjstart, totjend, reljstart, reljend


      call mppsync()
! ----------------------------------------
      IsBackup =  (datemean.eq.dateTo)
      if (lwp) write(*,*) 'diadump IsBackup = ',IsBackup, ' group ' ,FREQ_GROUP
! ----------------------------------------
      bkpname  = 'ave.20111231-15:30:00.N1p.nc.bkp'
      if (IsBackup) then
         forcing_file   = 'AVE_PHYS/ave.'//datemean//'.phys.nc.bkp'
      else
         forcing_file   = 'AVE_PHYS/ave.'//datemean//'.phys.nc'
      endif


      SELECT CASE (FREQ_GROUP)
        CASE (1) ; ave_counter=ave_counter_1 ; DIR='AVE_FREQ_1/'
        CASE (2) ; ave_counter=ave_counter_2 ; DIR='AVE_FREQ_2/'
      END SELECT

      jn_high = 0

! ******************  DIAGNOSTIC OUTPUT   2D *******************
      DO jn = 1, jptra_dia_b_2d

          if (.not.is_time_to_save(jn,FREQ_GROUP,2)) CYCLE
          if (FREQ_GROUP.eq.1) jn_high = jn_high+1
      if (rank == 0) then
! ******* rank 0 sets indexes of tot matrix where to place its own part

          iPd    = nldi
          iPe    = nlei
          jPd    = nldj
          jPe    = nlej
          istart = nimpp
          jstart = njmpp
          irange    = iPe - iPd + 1
          jrange    = jPe - jPd + 1
          totistart = istart + iPd - 1 ; totiend   = totistart + irange - 1
          totjstart = jstart + jPd - 1 ; totjend   = totjstart + jrange - 1
          relistart = 1 + iPd - 1      ; reliend   = relistart + irange - 1
          reljstart = 1 + jPd - 1      ; reljend   = reljstart + jrange - 1

          if (FREQ_GROUP.eq.1) then
          totbenIO2d(totistart:totiend, totjstart:totjend) = ben_DIA_2d_IO_HIGH(relistart:reliend,reljstart:reljend,jn_high)
          else
          totbenIO2d(totistart:totiend, totjstart:totjend) = ben_DIA_2d_IO(     relistart:reliend,reljstart:reljend,jn)
          endif

             do idrank = 1, size-1
! **************  rank 0 is receiving from the others their buffer  ****

                call MPI_RECV(jpi_rec    , 1,                 mpi_integer, idrank, 32,mpi_comm_world, status, ierr) !* first info to know where idrank is working
                call MPI_RECV(jpj_rec    , 1,                 mpi_integer, idrank, 33,mpi_comm_world, status, ierr)
                call MPI_RECV(istart     , 1,                 mpi_integer, idrank, 34,mpi_comm_world, status, ierr)
                call MPI_RECV(jstart     , 1,                 mpi_integer, idrank, 35,mpi_comm_world, status, ierr)
                call MPI_RECV(iPe        , 1,                 mpi_integer, idrank, 36,mpi_comm_world, status, ierr)
                call MPI_RECV(jPe        , 1,                 mpi_integer, idrank, 37,mpi_comm_world, status, ierr)
                call MPI_RECV(iPd        , 1,                 mpi_integer, idrank, 38,mpi_comm_world, status, ierr)
                call MPI_RECV(jPd        , 1                 ,mpi_integer, idrank, 39,mpi_comm_world, status, ierr)

                call MPI_RECV(buffben2D,jpi_rec*jpj_rec      ,mpi_real8,idrank, 40,mpi_comm_world, status, ierr)

! ******* rank 0 sets indexes of tot matrix where to place buffers of idrank
                irange    = iPe - iPd + 1
                jrange    = jPe - jPd + 1
                totistart = istart + iPd - 1 ; totiend   = totistart + irange - 1
                totjstart = jstart + jPd - 1 ; totjend   = totjstart + jrange - 1
                relistart = 1 + iPd - 1      ; reliend   = relistart + irange - 1
                reljstart = 1 + jPd - 1      ; reljend   = reljstart + jrange - 1

                do jj =totjstart,totjend ! only 2d vars
                 do ji =totistart,totiend
                  ind = (ji-totistart+ relistart )+ (jj-totjstart+ reljstart -1)*jpi_rec
                 totbenIO2d (ji,jj)= buffben2D(ind)
                 enddo
                enddo

             enddo !idrank = 1, size-1

      ELSE ! ranks 1 --> size-1


      if (FREQ_GROUP.eq.2) then
            do jj =1 , jpj
             do ji =1 , jpi
               ind           = ji + jpi * (jj-1)
               buffben2D (ind)= ben_DIA_2d_IO(ji,jj,jn)
              enddo
             enddo
      else
            do jj =1 , jpj
             do ji =1 , jpi
               ind           = ji + jpi * (jj-1)
               buffben2D (ind)= ben_DIA_2d_IO_high(ji,jj,jn_high)
              enddo
             enddo
      endif

              call MPI_SEND(jpi  , 1,mpi_integer, 0, 32, mpi_comm_world,ierr)
              call MPI_SEND(jpj  , 1,mpi_integer, 0, 33, mpi_comm_world,ierr)
              call MPI_SEND(nimpp, 1,mpi_integer, 0, 34, mpi_comm_world,ierr)
              call MPI_SEND(njmpp, 1,mpi_integer, 0, 35, mpi_comm_world,ierr)
              call MPI_SEND(nlei , 1,mpi_integer, 0, 36, mpi_comm_world,ierr)
              call MPI_SEND(nlej , 1,mpi_integer, 0, 37, mpi_comm_world,ierr)
              call MPI_SEND(nldi , 1,mpi_integer, 0, 38, mpi_comm_world,ierr)
              call MPI_SEND(nldj , 1,mpi_integer, 0, 39, mpi_comm_world,ierr)

             call MPI_SEND(buffben2D, jpi*jpj   ,mpi_real8, 0, 40, mpi_comm_world,ierr)

      ENDIF

      if (rank == 0) then
              var        =  dibnm_2d(jn)
              bkpname     = DIR//'ave.'//datemean//'.'//trim(var)//'.nc.bkp'
              dia_file_nc = DIR//'ave.'//datemean//'.'//trim(var)//'.nc'

              if (IsBackup) then
                 write(*,*) "trcdia ave_counter --> ", bkpname, ave_counter
                 CALL WRITE_AVE_2d_BKP_ben(bkpname,datefrom, dateTo,totbenIO2d, ave_counter);
              else
                 d2f2d_ben = REAL(totbenIO2d(:,:),4)
                 CALL WRITE_AVE_2d_ben(dia_file_nc,datefrom,dateTo, d2f2d_ben);

              endif
      end if ! if(rank == 0)

         if (.not.IsBackup) then
             if (FREQ_GROUP.eq.2) then
                ben_DIA_2d_IO(:,:,jn) = 0.
              else
                ben_DIA_2d_IO_HIGH(:,:,jn_high) = 0.
              endif
          endif

      ENDDO ! on jn


      jn_high = 0
! ******************  3D DIAGNOSTIC OUTPUT   *******************
      DO jn =1 , jptra_dia_b

          if (.not.is_time_to_save(jn,FREQ_GROUP,3)) CYCLE
          if (FREQ_GROUP.eq.1) jn_high = jn_high+1

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

      if (FREQ_GROUP.eq.1) then
      totbenIO (totistart:totiend, totjstart:totjend,:) = ben_DIA_IO_HIGH(relistart:reliend,reljstart:reljend, :,jn_high)
      else
      totbenIO (totistart:totiend, totjstart:totjend,:) = ben_DIA_IO(relistart:reliend, reljstart:reljend, :,jn) ! diagnostic from reaction model
      endif
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

                call MPI_RECV(buffben_dia ,jpi_rec*jpj_rec*jpk_ben,mpi_real8,idrank, 30,mpi_comm_world, status, ierr)

! ******* rank 0 sets indexes of tot matrix where to place buffers of idrank
                irange    = iPe - iPd + 1
                jrange    = jPe - jPd + 1
                totistart = istart + iPd - 1 ; totiend   = totistart + irange - 1
                totjstart = jstart + jPd - 1 ; totjend   = totjstart + jrange - 1
                relistart = 1 + iPd - 1      ; reliend   = relistart + irange - 1
                reljstart = 1 + jPd - 1      ; reljend   = reljstart + jrange - 1


                do jk =1 , jpk_ben ! 3d vars
                     do jj =totjstart,totjend
                        do ji =totistart,totiend
                           ind = (ji-totistart+ relistart )+ (jj-totjstart+ reljstart-1)*jpi_rec+(jk-1)*jpj_rec*jpi_rec
                           totbenIO(ji,jj,jk)= buffben_dia (ind)
                        enddo
                     enddo
                enddo


             enddo !idrank = 1, size-1


          else  ! IF LABEL 1,  if(rank == 0)

      if (FREQ_GROUP.eq.2) then
            do jk =1, jpk_ben
             do jj =1 , jpj
              do ji =1 , jpi
                  ind         =  ji + jpi * (jj-1) + jpi * jpj *(jk-1)
                  buffben_dia(ind) = ben_DIA_IO(ji,jj, jk,jn)
              enddo
             enddo
            enddo
      else
            do jk =1, jpk
             do jj =1 , jpj
              do ji =1 , jpi
                  ind         =  ji + jpi * (jj-1) + jpi * jpj *(jk-1)
                  buffben_dia(ind) = ben_DIA_IO_HIGH(ji,jj, jk,jn_high)
              enddo
             enddo
            enddo
      endif



              call MPI_SEND(jpi  , 1,mpi_integer, 0, 22, mpi_comm_world,ierr)
              call MPI_SEND(jpj  , 1,mpi_integer, 0, 23, mpi_comm_world,ierr)
              call MPI_SEND(nimpp, 1,mpi_integer, 0, 24, mpi_comm_world,ierr)
              call MPI_SEND(njmpp, 1,mpi_integer, 0, 25, mpi_comm_world,ierr)
              call MPI_SEND(nlei , 1,mpi_integer, 0, 26, mpi_comm_world,ierr)
              call MPI_SEND(nlej , 1,mpi_integer, 0, 27, mpi_comm_world,ierr)
              call MPI_SEND(nldi , 1,mpi_integer, 0, 28, mpi_comm_world,ierr)
              call MPI_SEND(nldj , 1,mpi_integer, 0, 29, mpi_comm_world,ierr)

              call MPI_SEND(buffben_dia  , jpi*jpj*jpk_ben,mpi_real8, 0, 30, mpi_comm_world,ierr)



      endif ! IF LABEL 1, if(rank == 0)
!************* END COLLECTING DATA  *****************

! *********** START WRITING **************************

      if (rank == 0) then
              var        =  dibnm(jn)
              bkpname     = DIR//'ave.'//datemean//'.'//trim(var)//'.nc.bkp'
              dia_file_nc = DIR//'ave.'//datemean//'.'//trim(var)//'.nc'

              if (IsBackup) then
                 !write(*,*) "trcdia ave_counter --> ", bkpname, ave_counter
                 CALL WRITE_AVE_BKP_ben(bkpname,datefrom, dateTo,totbenIO(:,:,:),ave_counter);
              else
                 d2f3d_ben = REAL(totbenIO(:,:,:),4)
                 CALL WRITE_AVE_ben(dia_file_nc,datefrom,dateTo, d2f3d_ben);

              endif


      end if ! IF LABEL 4  if(rank == 0)
         if (.not.IsBackup) then
             if (FREQ_GROUP.eq.2) then
                ben_DIA_IO(:,:,:,jn) = 0.
              else
                ben_DIA_IO_HIGH(:,:,:,jn_high) = 0.
              endif
          endif
      enddo  ! loop in jn


      CONTAINS

      LOGICAL FUNCTION IS_TIME_TO_SAVE(jn,FREQ_GROUP,ndims)
      IMPLICIT NONE

      integer jn, FREQ_GROUP,ndims

      IF (FREQ_GROUP.eq.2) then
         IS_TIME_TO_SAVE = .true.
      ELSE
        IS_TIME_TO_SAVE = .false.
        IF (ndims==3) then
            IF (dibhf(jn).eq.1)   IS_TIME_TO_SAVE = .true.
        ELSE
           IF (dibhf_2d(jn).eq.1) IS_TIME_TO_SAVE = .true.
        ENDIF
      ENDIF
      END FUNCTION IS_TIME_TO_SAVE

      end SUBROUTINE diawri_benthic
