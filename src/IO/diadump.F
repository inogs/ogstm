      SUBROUTINE diadump(datemean,datefrom,dateTo,FREQ_GROUP)
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
      INTEGER, INTENT(IN) :: FREQ_GROUP
      INTEGER ji, jj, jk, jn, jn_high
      INTEGER ind

      CHARACTER(LEN=42) forcing_file
      CHARACTER(LEN=43) bkpname
      CHARACTER(LEN=11) DIR
      logical IsBackup
      integer ave_counter


      CHARACTER(LEN=39) dia_file_nc
      CHARACTER(LEN=3)  var

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

      jn_high = 0

      SELECT CASE (FREQ_GROUP)
        CASE (1) ; ave_counter=ave_counter_1 ; DIR='AVE_FREQ_1/'
        CASE (2) ; ave_counter=ave_counter_2 ; DIR='AVE_FREQ_2/'
      END SELECT


!      PHYSICS FIRST!!
      if ( freq_ave_phys.eq.FREQ_GROUP) then
      ! *************** START COLLECTING DATA *****************************
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
          totistart = istart + iPd - 1 ; totiend   = totistart + irange - 1
          totjstart = jstart + jPd - 1 ; totjend   = totjstart + jrange - 1
          relistart = 1 + iPd - 1      ; reliend   = relistart + irange - 1
          reljstart = 1 + jPd - 1      ; reljend   = reljstart + jrange - 1

          totsnIO  (totistart:totiend, totjstart:totjend,:) = snIO    (relistart:reliend, reljstart:reljend, :)
          tottnIO  (totistart:totiend, totjstart:totjend,:) = tnIO    (relistart:reliend, reljstart:reljend, :)
          totvatmIO(totistart:totiend, totjstart:totjend)   = vatmIO  (relistart:reliend, reljstart:reljend)
          totempIO (totistart:totiend, totjstart:totjend)   = empIO   (relistart:reliend, reljstart:reljend)
          totqsrIO (totistart:totiend, totjstart:totjend)   = qsrIO   (relistart:reliend, reljstart:reljend)
          totunIO  (totistart:totiend, totjstart:totjend,:) = unIO    (relistart:reliend, reljstart:reljend, :)
          totbblxIO(totistart:totiend, totjstart:totjend)   = bblxIO  (relistart:reliend, reljstart:reljend)
          totvnIO  (totistart:totiend, totjstart:totjend,:) = vnIO    (relistart:reliend, reljstart:reljend, :)
          totbblyIO(totistart:totiend, totjstart:totjend)   = bblyIO  (relistart:reliend, reljstart:reljend)
          totwnIO  (totistart:totiend, totjstart:totjend,:) = wnIO    (relistart:reliend, reljstart:reljend, :)
          totavtIO (totistart:totiend, totjstart:totjend,:) = avtIO   (relistart:reliend, reljstart:reljend, :)
!          totdiaIO (totistart:totiend, totjstart:totjend,:,:) = tra_ppIO(relistart:reliend, reljstart:reljend, :,:) ! diagnostic from reaction model

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


      call MPI_RECV(buffsn  ,jpi_rec*jpj_rec*jpk          ,mpi_real8,idrank, 11,mpi_comm_world, status, ierr)
      call MPI_RECV(bufftn  ,jpi_rec*jpj_rec*jpk          ,mpi_real8,idrank, 12,mpi_comm_world, status, ierr)
      call MPI_RECV(buffvatm,jpi_rec*jpj_rec              ,mpi_real8,idrank, 13,mpi_comm_world, status, ierr)
      call MPI_RECV(buffemp ,jpi_rec*jpj_rec              ,mpi_real8,idrank, 14,mpi_comm_world, status, ierr)
      call MPI_RECV(buffqsr ,jpi_rec*jpj_rec              ,mpi_real8,idrank, 15,mpi_comm_world, status, ierr)
      call MPI_RECV(buffun  ,jpi_rec*jpj_rec*jpk          ,mpi_real8,idrank, 16,mpi_comm_world, status, ierr)
      call MPI_RECV(buffbblx,jpi_rec*jpj_rec              ,mpi_real8,idrank, 17,mpi_comm_world, status, ierr)
      call MPI_RECV(buffvn  ,jpi_rec*jpj_rec*jpk          ,mpi_real8,idrank, 18,mpi_comm_world, status, ierr)
      call MPI_RECV(buffbbly,jpi_rec*jpj_rec              ,mpi_real8,idrank, 19,mpi_comm_world, status, ierr)
      call MPI_RECV(buffwn  ,jpi_rec*jpj_rec*jpk          ,mpi_real8,idrank, 20,mpi_comm_world, status, ierr)
      call MPI_RECV(buffavt ,jpi_rec*jpj_rec*jpk          ,mpi_real8,idrank, 21,mpi_comm_world, status, ierr)
!      call MPI_RECV(buffpp  ,jpi_rec*jpj_rec*jpk*jptra_dia,mpi_real8,idrank, 22,mpi_comm_world, status, ierr)

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
                     totsnIO (ji,jj,jk)= buffsn (ind)
                     tottnIO (ji,jj,jk)= bufftn (ind)
                     totunIO (ji,jj,jk)= buffun (ind)
                     totvnIO (ji,jj,jk)= buffvn (ind)
                     totwnIO (ji,jj,jk)= buffwn (ind)
                     totavtIO(ji,jj,jk)= buffavt(ind)
                 enddo
                enddo
               enddo

!               do jn=1, jptra_dia
!                  do jk =1 , jpk ! 3d vars
!                     do jj =totjstart,totjend
!                        do ji =totistart,totiend
!                           ind = (ji-totistart+ relistart )+ (jj-totjstart+ reljstart-1)*jpi_rec+(jk-1)*jpj_rec*jpi_rec
!     &                           +(jn-1)*jpk*jpj_rec*jpi_rec
!                           totdiaIO(ji,jj,jk,jn)= buffpp (ind)
!                        enddo
!                     enddo
!                  enddo
!               enddo


               do jj =totjstart,totjend ! and 2d vars
                 do ji =totistart,totiend
                  ind = (ji-totistart+ relistart )+ (jj-totjstart+ reljstart -1)*jpi_rec
                  totvatmIO (ji,jj)= buffvatm(ind)
                  totempIO  (ji,jj)= buffemp (ind)
                  totqsrIO  (ji,jj)= buffqsr (ind)
                  totbblxIO (ji,jj)= buffbblx(ind)
                  totbblyIO (ji,jj)= buffbbly(ind)
             enddo
            enddo

          enddo !idrank = 1, size-1


      else  ! IF LABEL 1,  if(rank == 0)


!           do jn =1 , jptra_dia
!            do jk =1, jpk
!             do jj =1 , jpj
!              do ji =1 , jpi
!                  ind         =  ji + jpi * (jj-1) + jpi * jpj *(jk-1) + jpi*jpj*jpk *(jn-1)
!                  buffpp(ind) = tra_ppIO(ji,jj, jk,jn)
!              enddo
!             enddo
!            enddo
!           enddo

           do jk =1 , jpk
            do jj =1 , jpj
             do ji =1 , jpi
                   ind         =  ji + jpi * (jj-1) + jpi * jpj *(jk-1)
                   buffsn (ind)= snIO (ji,jj,jk)
                   bufftn (ind)= tnIO (ji,jj,jk)
                   buffun (ind)= unIO (ji,jj,jk)
                   buffvn (ind)= vnIO (ji,jj,jk)
                   buffwn (ind)= wnIO (ji,jj,jk)
                   buffavt(ind)= avtIO(ji,jj,jk)
              enddo
             enddo
            enddo

            do jj =1 , jpj
             do ji =1 , jpi
               ind           = ji + jpi * (jj-1)
               buffvatm (ind)= vatmIO(ji,jj)
               buffemp  (ind)= empIO (ji,jj)
               buffqsr  (ind)= qsrIO (ji,jj)
               buffbblx (ind)= bblxIO(ji,jj)
               buffbbly (ind)= bblyIO(ji,jj)
              enddo
             enddo


              call MPI_SEND(jpi  , 1,mpi_integer, 0, 1, mpi_comm_world,ierr)
              call MPI_SEND(jpj  , 1,mpi_integer, 0, 2, mpi_comm_world,ierr)
              call MPI_SEND(nimpp, 1,mpi_integer, 0, 3, mpi_comm_world,ierr)
              call MPI_SEND(njmpp, 1,mpi_integer, 0, 4, mpi_comm_world,ierr)
              call MPI_SEND(nlei , 1,mpi_integer, 0, 5, mpi_comm_world,ierr)
              call MPI_SEND(nlej , 1,mpi_integer, 0, 6, mpi_comm_world,ierr)
              call MPI_SEND(nldi , 1,mpi_integer, 0, 7, mpi_comm_world,ierr)
              call MPI_SEND(nldj , 1,mpi_integer, 0, 8, mpi_comm_world,ierr)

           call MPI_SEND(buffsn  , jpi*jpj*jpk  ,mpi_real8, 0, 11, mpi_comm_world,ierr)
           call MPI_SEND(bufftn  , jpi*jpj*jpk  ,mpi_real8, 0, 12, mpi_comm_world,ierr)
           call MPI_SEND(buffvatm, jpi*jpj      ,mpi_real8, 0, 13, mpi_comm_world,ierr)
           call MPI_SEND(buffemp , jpi*jpj      ,mpi_real8, 0, 14, mpi_comm_world,ierr)
           call MPI_SEND(buffqsr , jpi*jpj      ,mpi_real8, 0, 15, mpi_comm_world,ierr)
           call MPI_SEND(buffun  , jpi*jpj*jpk  ,mpi_real8, 0, 16, mpi_comm_world,ierr)
           call MPI_SEND(buffbblx, jpi*jpj      ,mpi_real8, 0, 17, mpi_comm_world,ierr)
           call MPI_SEND(buffvn  , jpi*jpj*jpk  ,mpi_real8, 0, 18, mpi_comm_world,ierr)
           call MPI_SEND(buffbbly, jpi*jpj      ,mpi_real8, 0, 19, mpi_comm_world,ierr)
           call MPI_SEND(buffwn  , jpi*jpj*jpk  ,mpi_real8, 0, 20, mpi_comm_world,ierr)
           call MPI_SEND(buffavt , jpi*jpj*jpk  ,mpi_real8, 0, 21, mpi_comm_world,ierr)




      endif ! IF LABEL 1, if(rank == 0)
!************* END COLLECTING DATA  *****************

! *********** START WRITING **************************

      if(rank == 0) then ! IF LABEL 4,
         if (IsBackup) then

            call PhysDump_bkp(forcing_file, datefrom, dateTo,ave_counter)
          else
            call PhysDump(forcing_file, datefrom, dateTo)
         endif
      endif

      endif !if ( freq_ave_phys.eq.FREQ_GROUP)


! ******************  DIAGNOSTIC OUTPUT   *******************
      DO jn =1 , jptra_dia

          if (.not.is_time_to_save(jn,FREQ_GROUP)) CYCLE

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
             jn_high = jn_high+1
      tottrnIO (totistart:totiend, totjstart:totjend,:) = tra_ppIO_HIGH(relistart:reliend,reljstart:reljend, :,jn_high)
             else
             tottrnIO (totistart:totiend, totjstart:totjend,:) = tra_ppIO(relistart:reliend, reljstart:reljend, :,jn) ! diagnostic from reaction model
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

      if (FREQ_GROUP.eq.2) then
            do jk =1, jpk
             do jj =1 , jpj
              do ji =1 , jpi
                  ind         =  ji + jpi * (jj-1) + jpi * jpj *(jk-1)
                  buffpp(ind) = tra_ppIO(ji,jj, jk,jn)
              enddo
             enddo
            enddo
      else
            do jk =1, jpk
             do jj =1 , jpj
              do ji =1 , jpi
                  ind         =  ji + jpi * (jj-1) + jpi * jpj *(jk-1)
                  buffpp(ind) = tra_ppIO_HIGH(ji,jj, jk,jn_high)
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

              call MPI_SEND(buffpp  , jpi*jpj*jpk,mpi_real8, 0, 30, mpi_comm_world,ierr)



      endif ! IF LABEL 1, if(rank == 0)
!************* END COLLECTING DATA  *****************

! *********** START WRITING **************************

      if (rank == 0) then
              var        =  dia(jn)
              bkpname     = DIR//'ave.'//datemean//'.'//var//'.nc.bkp'
              dia_file_nc = DIR//'ave.'//datemean//'.'//var//'.nc'

              if (IsBackup) then
                 !write(*,*) "trcdia ave_counter --> ", bkpname, ave_counter
                 CALL WRITE_AVE_BKP(bkpname,datefrom, dateTo,tottrnIO(:,:,:),ave_counter);
              else
                 d2f3d = REAL(tottrnIO(:,:,:),4)
                 CALL WRITE_AVE(dia_file_nc,datefrom,dateTo, d2f3d);

              endif


      end if ! IF LABEL 4  if(rank == 0)
         if (.not.IsBackup) then
             if (FREQ_GROUP.eq.2) then
                tra_ppIO(:,:,:,jn) = 0.
              else
                tra_ppIO_HIGH(:,:,:,jn_high) = 0.
              endif
          endif
      enddo  ! loop in jn

      if ((.not.IsBackup).and.( freq_ave_phys.eq.FREQ_GROUP) ) then

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
      endif


      CONTAINS

      LOGICAL FUNCTION IS_TIME_TO_SAVE(jn,FREQ_GROUP)
      IMPLICIT NONE

      integer jn, FREQ_GROUP

      IF (FREQ_GROUP.eq.2) then
         IS_TIME_TO_SAVE = .true.
      ELSE
        IF (highfreq_tra_dia(jn).eq.1) then
            IS_TIME_TO_SAVE = .true.
        ELSE
           IS_TIME_TO_SAVE = .false.
         ENDIF
      ENDIF
      END FUNCTION IS_TIME_TO_SAVE

      end SUBROUTINE diadump
