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
      use mpi
      USE ogstm_mpi_module


      IMPLICIT NONE

! local declarations
! ==================
      CHARACTER(LEN=17), INTENT(IN) :: datemean, dateFrom, dateTo
      INTEGER, INTENT(IN) :: FREQ_GROUP ! 1 = HIGH FREQ, 2 = LOW FREQ
      double precision ::  Miss_val =1.e20

      INTEGER jk,jj,ji, jn, jn_high,i,j,k
      INTEGER ind, i_contribution, j_contribution
      double precision :: elapsed_time

      CHARACTER(LEN=56) output_file_nc  ! AVE_FREQ_1/ave.20091231-12:00:00.P1n.nc
      CHARACTER(LEN=20)  var
      CHARACTER(LEN=60) bkpname
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
        CASE (1) 
       elapsed_time=elapsed_time_1
       DIR='AVE_FREQ_1/'
        CASE (2) 
       elapsed_time=elapsed_time_2
       DIR='AVE_FREQ_2/'
      END SELECT
      if (lwp) tottrnIO = Miss_val
!  Ghost Shells - Manual s Indexes-

      DO jn=1,jptra ! DO LABEL 5

      if (.not.is_time_to_save(jn,FREQ_GROUP)) CYCLE
      if (FREQ_GROUP.eq.1) jn_high = jn_high+1

      var    = ctrcnm(jn)
      output_file_nc = DIR//'ave.'//datemean//'.'//trim(var)//'.nc'
      bkpname        = DIR//'ave.'//datemean//'.'//trim(var)//'.nc.bkp'

! *************** START COLLECTING DATA *****************************
      if(myrank == 0) then                    ! IF LABEL 1

! ******* myrank 0 sets indexes of tot matrix where to place its own part
          iPd    = nldi
          iPe    = nlei
          jPd    = nldj
          jPe    = nlej
          istart = nimpp
          jstart = njmpp
          irange    = iPe - iPd + 1
          jrange    = jPe - jPd + 1
          totistart = istart + iPd - 1
          totiend  = totistart + irange - 1
          totjstart = jstart + jPd - 1
          totjend  = totjstart + jrange - 1
          relistart = 1 + iPd - 1     
          reliend  = relistart + irange - 1
          reljstart = 1 + jPd - 1     
          reljend  = reljstart + jrange - 1


!***** START ASSEMBLING ***  myrank 0 puts its tracer part in the tot matrix
          if (FREQ_GROUP.eq.1) then
          tottrnIO(:,totjstart:totjend,totistart:totiend)= traIO_HIGH(:,reljstart:reljend,relistart:reliend,jn_high)
          else
          tottrnIO(:,totjstart:totjend,totistart:totiend)= traIO (:,reljstart:reljend,relistart:reliend,jn)
          endif


           do idrank = 1,mpi_glcomm_size-1
! **************  myrank 0 is receiving from the others their buffer  ****
           call MPI_RECV(jpi_rec    , 1,                 mpi_integer, idrank, 1,mpi_comm_world, status, ierr) !* first info to know where idrank is working
           call MPI_RECV(jpj_rec    , 1,                 mpi_integer, idrank, 2,mpi_comm_world, status, ierr)
           call MPI_RECV(istart     , 1,                 mpi_integer, idrank, 3,mpi_comm_world, status, ierr)
           call MPI_RECV(jstart     , 1,                 mpi_integer, idrank, 4,mpi_comm_world, status, ierr)
           call MPI_RECV(iPe        , 1,                 mpi_integer, idrank, 5,mpi_comm_world, status, ierr)
           call MPI_RECV(jPe        , 1,                 mpi_integer, idrank, 6,mpi_comm_world, status, ierr)
           call MPI_RECV(iPd        , 1,                 mpi_integer, idrank, 7,mpi_comm_world, status, ierr)
           call MPI_RECV(jPd        , 1                 ,mpi_integer, idrank, 8,mpi_comm_world, status, ierr)
           call MPI_RECV(bufftrn    ,jpi_rec*jpj_rec*jpk,  mpi_real8, idrank, 9,mpi_comm_world, status, ierr) ! ** then tracer buffer
! ******* myrank 0 sets indexes of tot matrix where to place buffers of idrank
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
! **** ASSEMBLING *** myrank 0 puts in tot matrix buffer received by idrank
          do ji =totistart,totiend
            i_contribution   = jpk*jpj_rec*(ji-1-totistart+ relistart)
            do jj =totjstart,totjend
              j_contribution = jpk*(jj-1-totjstart+ reljstart)
              do jk =1, jpk
                  ind = jk + j_contribution + i_contribution
                  tottrnIO(jk,jj,ji)= bufftrn(ind)
              enddo
            enddo
          enddo
!      ciclo vecchia maniera
!          do jk =1 , jpk
!           do jj =totjstart,totjend
!             do ji =totistart,totiend
!                 ind = (ji-totistart+ relistart)+ (jj-totjstart+ reljstart-1)*jpi_rec+(jk-1)*jpj_rec*jpi_rec
!                 tottrnIO(jk,jj,ji)= bufftrn   (ind)
!             enddo
!            enddo
!           enddo
       enddo !idrank = 1, size-1
   else  ! IF LABEL 1,  if(myrank == 0)
! **** work of the other ranks
! ****** 1. load  inf buffer their IO matrices
     if (FREQ_GROUP.eq.2) then
        do ji =1 , jpi
          i_contribution= jpk*jpj * (ji - 1 )
         do jj =1 , jpj
          j_contribution=jpk*(jj-1)
          do jk =1 , jpk
            ind            =  jk + j_contribution + i_contribution
            bufftrn   (ind)= traIO( jk,jj,ji,jn)
          enddo
         enddo
       enddo
     else ! FREQ_GROUP.eq.1
        do ji =1 , jpi
          i_contribution= jpk*jpj * (ji - 1 )
         do jj =1 , jpj
          j_contribution=jpk*(jj-1)
          do jk =1 , jpk
            ind            =  jk + j_contribution + i_contribution
            bufftrn   (ind)= traIO_HIGH( jk,jj,ji,jn_high)
          enddo
         enddo
       enddo
     endif
! ******  2.send buffer to myrank 0
       call MPI_SEND(jpi  , 1,mpi_integer, 0, 1, mpi_comm_world,ierr)
       call MPI_SEND(jpj  , 1,mpi_integer, 0, 2, mpi_comm_world,ierr)
       call MPI_SEND(nimpp, 1,mpi_integer, 0, 3, mpi_comm_world,ierr)
       call MPI_SEND(njmpp, 1,mpi_integer, 0, 4, mpi_comm_world,ierr)
       call MPI_SEND(nlei , 1,mpi_integer, 0, 5, mpi_comm_world,ierr)
       call MPI_SEND(nlej , 1,mpi_integer, 0, 6, mpi_comm_world,ierr)
       call MPI_SEND(nldi , 1,mpi_integer, 0, 7, mpi_comm_world,ierr)
       call MPI_SEND(nldj , 1,mpi_integer, 0, 8, mpi_comm_world,ierr)
       call MPI_SEND(bufftrn, jpi * jpj * jpk,MPI_DOUBLE_PRECISION, 0, 9, mpi_comm_world,ierr)
   endif ! IF LABEL 1, if(myrank == 0)

!************* END COLLECTING DATA  *****************



! *********** START WRITING **************************
      if(myrank == 0) then ! IF LABEL 4

        if (IsBackup) then
          CALL WRITE_AVE_BKP(bkpname,var,datefrom, dateTo,tottrnIO,elapsed_time, deflate_ave, deflate_level_ave)

        else

          CALL WRITE_AVE(output_file_nc,var,datefrom, dateTo, tottrnIO, deflate_ave, deflate_level_ave)
      
         endif


      end if ! IF LABEL 4  if(myrank == 0)

! 200     FORMAT(' ',A3,I4,I4,I4,I4,D30.23)



      if (.not.IsBackup) then
         if (FREQ_GROUP.eq.2) then
            traIO(:,:,:,jn) = 0.  !      we reset matrix for new average
         else
            traIO_HIGH(:,:,:,jn_high) = 0.
         endif
      endif

      END DO ! do jn=1,jptra , DO LABEL 5

      ! stop




      if(myrank == 0) WRITE(numout,*) '**** trcdit : write NetCDF passive tracer concentration'

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
