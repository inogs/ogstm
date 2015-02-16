      SUBROUTINE trcdit(datemean,datefrom,dateTo,FREQ_GROUP)
!---------------------------------------------------------------------
!
!                       ROUTINE trcdit
!                     ******************
!
!  Purpose :
!  ---------
!     Standard output of passive tracer : concentration fields



      USE calendar
      USE myalloc
      USE IO_mem
      USE FN_mem
      USE TIME_MANAGER
#ifdef key_mpp
      USE myalloc_mpp
#endif

      IMPLICIT NONE

! local declarations
! ==================
      CHARACTER(LEN=17), INTENT(IN) :: datemean, dateFrom, dateTo
      INTEGER, INTENT(IN) :: FREQ_GROUP ! 1 = HIGH FREQ, 2 = LOW FREQ


      INTEGER ji, jj, jk, jn, jn_high
      INTEGER ind
      INTEGER ave_counter

      CHARACTER(LEN=39) output_file_nc  ! AVE_FREQ_1/ave.20091231-12:00:00.P1n.nc
      CHARACTER(LEN=3)  var
      CHARACTER(LEN=43) bkpname
      CHARACTER(LEN=11) DIR
      logical IsBackup

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
      call mppsync()

      jn_high = 0

      SELECT CASE (FREQ_GROUP)
        CASE (1) ; ave_counter=ave_counter_1 ; DIR='AVE_FREQ_1/'
        CASE (2) ; ave_counter=ave_counter_2 ; DIR='AVE_FREQ_2/'
      END SELECT

!  Ghost Shells - Manual s Indexes-

      DO jn=1,jptra ! DO LABEL 5

      if (.not.is_time_to_save(jn,FREQ_GROUP)) CYCLE
      if (FREQ_GROUP.eq.1) jn_high = jn_high+1

      var    = ctrcnm(jn)
      output_file_nc = DIR//'ave.'//datemean//'.'//var//'.nc'
      bkpname        = DIR//'ave.'//datemean//'.'//var//'.nc.bkp'

! *************** START COLLECTING DATA *****************************
      if(rank == 0) then                    ! IF LABEL 1

! ******* rank 0 sets indexes of tot matrix where to place its own part
          iPd    = nldi
          iPe    = nlei
          jPd    = nldj
          jPe    = nlej
          istart = nimpp
          jstart = njmpp
          irange    = iPe - iPd + 1
          jrange    = jPe - jPd + 1
          totistart = istart + iPd - 1; totiend  = totistart + irange - 1
          totjstart = jstart + jPd - 1; totjend  = totjstart + jrange - 1
          relistart = 1 + iPd - 1     ; reliend  = relistart + irange - 1
          reljstart = 1 + jPd - 1     ; reljend  = reljstart + jrange - 1


!***** START ASSEMBLING ***  rank 0 puts its tracer part in the tot matrix
          if (FREQ_GROUP.eq.1) then
          tottrnIO(totistart:totiend, totjstart:totjend,:)= traIO_HIGH(relistart:reliend, reljstart:reljend,:,jn_high)
          else
          tottrnIO(totistart:totiend, totjstart:totjend,:)= traIO (relistart:reliend, reljstart:reljend,:,jn)
          endif


          do idrank = 1, size-1
! **************  rank 0 is receiving from the others their buffer  ****
              call MPI_RECV(jpi_rec    , 1,                 mpi_integer, idrank, 1,mpi_comm_world, status, ierr) !* first info to know where idrank is working
              call MPI_RECV(jpj_rec    , 1,                 mpi_integer, idrank, 2,mpi_comm_world, status, ierr)
              call MPI_RECV(istart     , 1,                 mpi_integer, idrank, 3,mpi_comm_world, status, ierr)
              call MPI_RECV(jstart     , 1,                 mpi_integer, idrank, 4,mpi_comm_world, status, ierr)
              call MPI_RECV(iPe        , 1,                 mpi_integer, idrank, 5,mpi_comm_world, status, ierr)
              call MPI_RECV(jPe        , 1,                 mpi_integer, idrank, 6,mpi_comm_world, status, ierr)
              call MPI_RECV(iPd        , 1,                 mpi_integer, idrank, 7,mpi_comm_world, status, ierr)
              call MPI_RECV(jPd        , 1                 ,mpi_integer, idrank, 8,mpi_comm_world, status, ierr)
              call MPI_RECV(bufftrn    ,jpi_rec*jpj_rec*jpk,  mpi_real8, idrank, 9,mpi_comm_world, status, ierr) ! ** then tracer buffer

! ******* rank 0 sets indexes of tot matrix where to place buffers of idrank
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
! **** ASSEMBLING *** rank 0 puts in tot matrix buffer received by idrank
              do jk =1 , jpk
               do jj =totjstart,totjend
                 do ji =totistart,totiend
                     ind = (ji-totistart+ relistart)+ (jj-totjstart+ reljstart-1)*jpi_rec+(jk-1)*jpj_rec*jpi_rec
                     tottrnIO(ji,jj,jk)= bufftrn   (ind)
                 enddo
                enddo
               enddo

          enddo !idrank = 1, size-1


      else  ! IF LABEL 1,  if(rank == 0)
! **** work of the other ranks
! ****** 1. load  inf buffer their IO matrices

        if (FREQ_GROUP.eq.2) then
           do jk =1 , jpk
            do jj =1 , jpj
             do ji =1 , jpi
               ind            =  ji + jpi * (jj-1) + jpi * jpj *(jk-1)
               bufftrn   (ind)= traIO( ji,jj,jk,jn)
             enddo
            enddo
          enddo
        else ! FREQ_GROUP.eq.1
           do jk =1 , jpk
            do jj =1 , jpj
             do ji =1 , jpi
               ind            =  ji + jpi * (jj-1) + jpi * jpj *(jk-1)
               bufftrn   (ind)= traIO_HIGH( ji,jj,jk,jn_high)
             enddo
            enddo
          enddo
        endif




! ******  2.send buffer to rank 0
          call MPI_SEND(jpi  , 1,mpi_integer, 0, 1, mpi_comm_world,ierr)
          call MPI_SEND(jpj  , 1,mpi_integer, 0, 2, mpi_comm_world,ierr)
          call MPI_SEND(nimpp, 1,mpi_integer, 0, 3, mpi_comm_world,ierr)
          call MPI_SEND(njmpp, 1,mpi_integer, 0, 4, mpi_comm_world,ierr)
          call MPI_SEND(nlei , 1,mpi_integer, 0, 5, mpi_comm_world,ierr)
          call MPI_SEND(nlej , 1,mpi_integer, 0, 6, mpi_comm_world,ierr)
          call MPI_SEND(nldi , 1,mpi_integer, 0, 7, mpi_comm_world,ierr)
          call MPI_SEND(nldj , 1,mpi_integer, 0, 8, mpi_comm_world,ierr)
          call MPI_SEND(bufftrn, jpi * jpj * jpk,mpi_real8  , 0, 9, mpi_comm_world,ierr)



      endif ! IF LABEL 1, if(rank == 0)

!************* END COLLECTING DATA  *****************



! *********** START WRITING **************************
      if(rank == 0) then ! IF LABEL 4


        if (IsBackup) then
          CALL WRITE_AVE_BKP(bkpname,datefrom, dateTo,tottrnIO,ave_counter);

        else
          d2f3d = REAL(tottrnIO,4)
          CALL WRITE_AVE(output_file_nc,datefrom, dateTo, d2f3d);
         endif


      end if ! IF LABEL 4  if(rank == 0)




      if (.not.IsBackup) then
         if (FREQ_GROUP.eq.2) then
            traIO(:,:,:,jn) = 0.  !      we reset matrix for new average
         else
            traIO_HIGH(:,:,:,jn_high) = 0.
         endif
      endif

      END DO ! do jn=1,jptra , DO LABEL 5





      if(rank == 0) WRITE(numout,*) '**** trcdit : write NetCDF passive tracer concentration'

      CONTAINS

      LOGICAL FUNCTION IS_TIME_TO_SAVE(jn,FREQ_GROUP)
      IMPLICIT NONE

      integer jn, FREQ_GROUP

      IF (FREQ_GROUP.eq.2) then
         IS_TIME_TO_SAVE = .true.
      ELSE
        IF (ctr_hf(jn).eq.1) then
            IS_TIME_TO_SAVE = .true.
        ELSE
           IS_TIME_TO_SAVE = .false.
         ENDIF
      ENDIF
      END FUNCTION IS_TIME_TO_SAVE

      END SUBROUTINE trcdit





