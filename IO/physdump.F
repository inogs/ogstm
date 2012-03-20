      SUBROUTINE physdump(datemean,datefrom,dateTo)
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

      REAL(4) lon(jpiglo), lat(jpjglo), dep(jpk)

      INTEGER ji, jj, jk, jn ,jn_nudg, jl
      INTEGER ind, AC

      CHARACTER*32 forcing2d_file, forcing3d_file
      CHARACTER*33 bkpname
      logical IsBackup


      CHARACTER*29 dia_file_nc
      CHARACTER*4  var ! , var_nc

      INTEGER :: ir, ir2d, ir3d
      INTEGER idrank, ierr, istart, jstart, iPe, iPd, jPe, jPd, status(MPI_STATUS_SIZE)
      INTEGER irange, jrange
      INTEGER totistart, totiend, relistart, reliend
      INTEGER totjstart, totjend, reljstart, reljend


      call mppsync()
! ----------------------------------------
      IsBackup =  (datemean.eq.dateTo)
      if (lwp) write(*,*) 'physdump IsBackup = ',IsBackup
! ----------------------------------------
      bkpname  = 'ave.20111231-15:30:00.N1p.nc.bkp'//CHAR(0)
      AC       = ave_counter


      forcing3D_file = 'ave.'//datemean//'.grid3D.dat'
      forcing2D_file = 'ave.'//datemean//'.grid2D.dat'

      var   (4:4)    = CHAR(0)
!      var_nc(4:4)    = CHAR(0)


      ! *************** START COLLECTING DATA *****************************
      if (rank == 0) then                    ! IF LABEL 1
      ir2d=jpiglo*jpjglo*4
      ir3d=jpiglo*jpjglo*jpk*4

! ******* rank 0 sets indexes of tot matrix where to place its own part

          iPd    = nldi
          iPe    = nlei
          jPd    = nldj
          jPe    = nlej
          istart = nimpp
          jstart = njmpp
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



          write(*,*) forcing3D_file
          open(4446, file=forcing3D_file, status='unknown',form='unformatted',access='direct', recl=ir3d) ! Output of mean Forcing Fields grid.T,.U,.V,.W
          open(4447, file=forcing2D_file, status='unknown',form='unformatted',access='direct', recl=ir2d)

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
          totdiaIO (totistart:totiend, totjstart:totjend,:,:) = tra_ppIO(relistart:reliend, reljstart:reljend, :,:) ! diagnostic from reaction model

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
      call MPI_RECV(buffpp  ,jpi_rec*jpj_rec*jpk*jptra_dia,mpi_real8,idrank, 22,mpi_comm_world, status, ierr)

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

               do jn=1, jptra_dia
                  do jk =1 , jpk ! 3d vars
                     do jj =totjstart,totjend
                        do ji =totistart,totiend
                           ind = (ji-totistart+ relistart )+ (jj-totjstart+ reljstart-1)*jpi_rec+(jk-1)*jpj_rec*jpi_rec
     &                           +(jn-1)*jpk*jpj_rec*jpi_rec
                           totdiaIO(ji,jj,jk,jn)= buffpp (ind)
                        enddo
                     enddo
                  enddo
               enddo


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


           do jn =1 , jptra_dia
            do jk =1, jpk
             do jj =1 , jpj
              do ji =1 , jpi
                  ind         =  ji + jpi * (jj-1) + jpi * jpj *(jk-1) + jpi*jpj*jpk *(jn-1)
                  buffpp(ind) = tra_ppIO(ji,jj, jk,jn)
              enddo
             enddo
            enddo
           enddo

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
           call MPI_SEND(buffpp  , jpi*jpj*jpk*jptra_dia,mpi_real8, 0, 22, mpi_comm_world,ierr)



      endif ! IF LABEL 1, if(rank == 0)


      if(rank == 0) then ! IF LABEL 4,


            DO ji=1, jpiglo ; lon(ji)= REAL(glamt(1,1),4)+(ji-1)*0.125 ; ENDDO
            DO jj=1, jpjglo ; lat(jj)= REAL(gphit(1,1),4)+(jj-1)*0.125 ; ENDDO
            DO jk=1, jpk    ; dep(jk)= REAL(gdept(jk),4)               ; ENDDO


            write(4446,rec=1) (((REAL(totsnIO (ji,jj,jk),4), ji=1,jpiglo), jj=1,jpjglo), jk=1,jpk)
            write(4446,rec=2) (((REAL(tottnIO (ji,jj,jk),4), ji=1,jpiglo), jj=1,jpjglo), jk=1,jpk)
            write(4446,rec=3) (((REAL(totunIO (ji,jj,jk),4), ji=1,jpiglo), jj=1,jpjglo), jk=1,jpk)
            write(4446,rec=4) (((REAL(totvnIO (ji,jj,jk),4), ji=1,jpiglo), jj=1,jpjglo), jk=1,jpk)
            write(4446,rec=5) (((REAL(totwnIO (ji,jj,jk),4), ji=1,jpiglo), jj=1,jpjglo), jk=1,jpk)
            write(4446,rec=6) (((REAL(totavtIO(ji,jj,jk),4), ji=1,jpiglo), jj=1,jpjglo), jk=1,jpk)
            write(4447,rec=1) ((REAL(totvatmIO(ji,jj)   ,4), ji=1,jpiglo), jj=1,jpjglo)
            write(4447,rec=2) ((REAL(totempIO (ji,jj)   ,4), ji=1,jpiglo), jj=1,jpjglo)
            write(4447,rec=3) ((REAL(totqsrIO (ji,jj)   ,4), ji=1,jpiglo), jj=1,jpjglo)
            write(4447,rec=4) ((REAL(totbblxIO(ji,jj)   ,4), ji=1,jpiglo), jj=1,jpjglo)
            write(4447,rec=5) ((REAL(totbblyIO(ji,jj)   ,4), ji=1,jpiglo), jj=1,jpjglo)


! ******************  DIAGNOSTIC OUTPUT   *******************
           do jn =1 , jptra_dia
              var(1:3) = dia(jn)
              bkpname    = 'ave.'//datemean//'.'//dia(jn)//'.nc.bkp'//CHAR(0)
              dia_file_nc ='ave.'//datemean//'.'//dia(jn)//'.nc'//CHAR(0)

              if (IsBackup) then
                 CALL IOOGSNC_AVE_BKP(bkpname,var,jpiglo,jpjglo,jpk,lon,lat,dep,datefrom,dateTo,totdiaIO(:,:,:,jn),AC)
              else
                 d2f3d= REAL(totdiaIO(:,:,:,jn),4)     !   conversion from double to float
                 CALL IOOGSNC_AVE_3D2(dia_file_nc, var,jpiglo,jpjglo,jpk,lon,lat,dep,datefrom, dateTo,d2f3d)
              endif
           enddo


!           d2f3d= REAL(totdiaIO(:,:,:,1),4)     !   conversion from double to float
!           var(1:3) = 'ppg'
!           CALL IOOGSNC_AVE_3D2(dia1_file_nc, var,jpiglo,jpjglo,jpk,lon,lat,dep,datefrom, dateTo,d2f3d)
!
! *****************  output NET  PRIMARY  PRODUCTION ********************
!
!           d2f3d= REAL(totdiaIO(:,:,:,2),4)    !   conversion from double to float
!           var(1:3) = 'ppn'
!           CALL IOOGSNC_AVE_3D2(dia2_file_nc, var,jpiglo,jpjglo,jpk,lon,lat,dep,datefrom, dateTo,d2f3d)
!


            close(4446)
            close(4447)


      end if ! IF LABEL 4  if(rank == 0)



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





      end SUBROUTINE physdump
