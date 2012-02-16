

      SUBROUTINE trcdit(datemean,datefrom,dateTo)
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

      REAL(4) lon(jpiglo), lat(jpjglo), dep(jpk)

      INTEGER ji, jj, jk, jn,jn_nudg, jl
      INTEGER ind

!      CHARACTER*32 forcing2d_file, forcing3d_file
      CHARACTER*29  output_file_FN

      CHARACTER*29 output_file_nc !,ppg_file_nc, ppn_file_nc,ppb_file_nc  ! ave.20091231-12:00:00.P1n.nc
      CHARACTER*4  var, var_nc


!----------------------------------------------------------------------
! statement functions
! ===================

#include "stafun.h"

      INTEGER :: ir, ir2d, ir3d


      INTEGER idrank, ierr, istart, jstart, iPe, iPd, jPe, jPd, status(MPI_STATUS_SIZE)
      INTEGER irange, jrange
      INTEGER totistart, totiend, relistart, reliend
      INTEGER totjstart, totjend, reljstart, reljend


      call mppsync()
      ave_counter = 0   !  we reset the counter


!      forcing3D_file = 'ave.'//datemean//'.grid3D.dat'
!      forcing2D_file = 'ave.'//datemean//'.grid2D.dat'
!
!      ppg_file_nc    = 'ave.'//datemean//'.ppg.nc'//CHAR(0)
!      ppn_file_nc    = 'ave.'//datemean//'.ppn.nc'//CHAR(0)
!      ppb_file_nc    = 'ave.'//datemean//'.ppb.nc'//CHAR(0)

      var   (4:4)    = CHAR(0)
      var_nc(4:4)    = CHAR(0)




!  Ghost Shells - Manual s Indexes-

      DO jn=1,jptra ! DO LABEL 5

      var_nc(1:3)    = ctrcnm(jn)(1:3)
      output_file_nc = 'ave.'//datemean//'.'//var_nc(1:3)//'.nc'//CHAR(0)
      output_file_FN = 'fne.'//datemean//'.'//var_nc(1:3)//'.dat'

! *************** START COLLECTING DATA *****************************
      if(rank == 0) then                    ! IF LABEL 1
          ir  =jpiglo*jpjglo*jpk*4
          ir3d=jpiglo*jpjglo*jpk*4
          ir2d=jpiglo*jpjglo*4

          open(4445, file=output_file_FN, status='unknown',form='unformatted',access='direct', recl=ir)

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


!***** START ASSEMBLING ***  rank 0 puts its tracer part in the tot matrix
          tottrnIO(totistart:totiend, totjstart:totjend,:)= traIO (relistart:reliend, reljstart:reljend,:,jn)
          tottrnFN(totistart:totiend, totjstart:totjend,:)= tra_FN(relistart:reliend, reljstart:reljend,:,jn)

!          IF (jn .EQ. 1) THEN ! the same for PHYSICAL PART once and for all
!              open(4446, file=forcing3D_file, status='unknown',form='unformatted',access='direct', recl=ir3d) ! Output of mean Forcing Fields grid.T,.U,.V,.W
!              open(4447, file=forcing2D_file, status='unknown',form='unformatted',access='direct', recl=ir2d)
!
!              totsnIO  (totistart:totiend, totjstart:totjend,:) = snIO    (relistart:reliend, reljstart:reljend, :)
!              tottnIO  (totistart:totiend, totjstart:totjend,:) = tnIO    (relistart:reliend, reljstart:reljend, :)
!              totvatmIO(totistart:totiend, totjstart:totjend)   = vatmIO  (relistart:reliend, reljstart:reljend)
!              totempIO (totistart:totiend, totjstart:totjend)   = empIO   (relistart:reliend, reljstart:reljend)
!              totqsrIO (totistart:totiend, totjstart:totjend)   = qsrIO   (relistart:reliend, reljstart:reljend)
!              totunIO  (totistart:totiend, totjstart:totjend,:) = unIO    (relistart:reliend, reljstart:reljend, :)
!              totbblxIO(totistart:totiend, totjstart:totjend)   = bblxIO  (relistart:reliend, reljstart:reljend)
!              totvnIO  (totistart:totiend, totjstart:totjend,:) = vnIO    (relistart:reliend, reljstart:reljend, :)
!              totbblyIO(totistart:totiend, totjstart:totjend)   = bblyIO  (relistart:reliend, reljstart:reljend)
!              totwnIO  (totistart:totiend, totjstart:totjend,:) = wnIO    (relistart:reliend, reljstart:reljend, :)
!              totavtIO (totistart:totiend, totjstart:totjend,:) = avtIO   (relistart:reliend, reljstart:reljend, :)
!              totppgIO (totistart:totiend, totjstart:totjend,:) = tra_ppIO(relistart:reliend, reljstart:reljend, :,1)
!              totppnIO (totistart:totiend, totjstart:totjend,:) = tra_ppIO(relistart:reliend, reljstart:reljend, :,2)
!              totppbIO (totistart:totiend, totjstart:totjend,:) = tra_ppIO(relistart:reliend, reljstart:reljend, :,3)
!          ENDIF


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
              call MPI_RECV(bufftrn_FN, jpi_rec*jpj_rec*jpk,  mpi_real8, idrank,10,mpi_comm_world, status, ierr)
!              IF (jn .EQ. 1) THEN  ! ********* receiving of physical data, once and for all
!                  call MPI_RECV(buffsn  ,jpi_rec*jpj_rec*jpk  ,mpi_real8,idrank, 11,mpi_comm_world, status, ierr)
!                  call MPI_RECV(bufftn  ,jpi_rec*jpj_rec*jpk  ,mpi_real8,idrank, 12,mpi_comm_world, status, ierr)
!                  call MPI_RECV(buffvatm,jpi_rec*jpj_rec      ,mpi_real8,idrank, 13,mpi_comm_world, status, ierr)
!                  call MPI_RECV(buffemp ,jpi_rec*jpj_rec      ,mpi_real8,idrank, 14,mpi_comm_world, status, ierr)
!                  call MPI_RECV(buffqsr ,jpi_rec*jpj_rec      ,mpi_real8,idrank, 15,mpi_comm_world, status, ierr)
!                  call MPI_RECV(buffun  ,jpi_rec*jpj_rec*jpk  ,mpi_real8,idrank, 16,mpi_comm_world, status, ierr)
!                  call MPI_RECV(buffbblx,jpi_rec*jpj_rec      ,mpi_real8,idrank, 17,mpi_comm_world, status, ierr)
!                  call MPI_RECV(buffvn  ,jpi_rec*jpj_rec*jpk  ,mpi_real8,idrank, 18,mpi_comm_world, status, ierr)
!                  call MPI_RECV(buffbbly,jpi_rec*jpj_rec      ,mpi_real8,idrank, 19,mpi_comm_world, status, ierr)
!                  call MPI_RECV(buffwn  ,jpi_rec*jpj_rec*jpk  ,mpi_real8,idrank, 20,mpi_comm_world, status, ierr)
!                  call MPI_RECV(buffavt ,jpi_rec*jpj_rec*jpk  ,mpi_real8,idrank, 21,mpi_comm_world, status, ierr)
!                  call MPI_RECV(buffpp  ,jpi_rec*jpj_rec*jpk*3,mpi_real8,idrank, 22,mpi_comm_world, status, ierr)
!             ENDIF
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
                     tottrnFN(ji,jj,jk)= bufftrn_FN(ind)
                 enddo
                enddo
               enddo

! ******* ASSEMBLING of PHYSICAL PART, once and for all
!              IF (jn .EQ.1 ) THEN
!
!                  do jk =1 , jpk ! 3d vars
!                   do jj =totjstart,totjend
!                     do ji =totistart,totiend
!                         ind = (ji-totistart+ relistart )+ (jj-totjstart+ reljstart-1)*jpi_rec+(jk-1)*jpj_rec*jpi_rec
!                         totsnIO (ji,jj,jk)= buffsn (ind)
!                         tottnIO (ji,jj,jk)= bufftn (ind)
!                         totunIO (ji,jj,jk)= buffun (ind)
!                         totvnIO (ji,jj,jk)= buffvn (ind)
!                         totwnIO (ji,jj,jk)= buffwn (ind)
!                         totavtIO(ji,jj,jk)= buffavt(ind)
!                         totppgIO(ji,jj,jk)= buffpp (ind)
!                         totppnIO(ji,jj,jk)= buffpp (ind+ 1*jpk*jpj_rec*jpi_rec)
!                         totppbIO(ji,jj,jk)= buffpp (ind+ 2*jpk*jpj_rec*jpi_rec)
!                     enddo
!                    enddo
!                   enddo
!
!                   do jj =totjstart,totjend ! and 2d vars
!                     do ji =totistart,totiend
!                      ind = (ji-totistart+ relistart )+ (jj-totjstart+ reljstart -1)*jpi_rec
!                      totvatmIO (ji,jj)= buffvatm(ind)
!                      totempIO  (ji,jj)= buffemp (ind)
!                      totqsrIO  (ji,jj)= buffqsr (ind)
!                      totbblxIO (ji,jj)= buffbblx(ind)
!                      totbblyIO (ji,jj)= buffbbly(ind)
!                 enddo
!                enddo
!
!              ENDIF


          enddo !idrank = 1, size-1


      else  ! IF LABEL 1,  if(rank == 0)
! **** work of the other ranks
! ****** 1. load  inf buffer their IO matrices


           do jk =1 , jpk
            do jj =1 , jpj
             do ji =1 , jpi
               ind            =  ji + jpi * (jj-1) + jpi * jpj *(jk-1)
               bufftrn   (ind)= traIO( ji,jj,jk,jn)
               bufftrn_FN(ind)= tra_FN(ji,jj,jk,jn)
             enddo
            enddo
          enddo


!        IF (jn .EQ. 1) THEN
!               do jn_nudg =1 , 3
!                do jk =1, jpk
!                 do jj =1 , jpj
!                  do ji =1 , jpi
!                      ind         =  ji + jpi * (jj-1) + jpi * jpj *(jk-1) + jpi*jpj*jpk *(jn_nudg-1)
!                      buffpp(ind) = tra_ppIO(ji,jj, jk,jn_nudg)
!                  enddo
!                 enddo
!                enddo
!               enddo
!
!               do jk =1 , jpk
!                do jj =1 , jpj
!                 do ji =1 , jpi
!                       ind         =  ji + jpi * (jj-1) + jpi * jpj *(jk-1)
!                       buffsn (ind)= snIO (ji,jj,jk)
!                       bufftn (ind)= tnIO (ji,jj,jk)
!                       buffun (ind)= unIO (ji,jj,jk)
!                       buffvn (ind)= vnIO (ji,jj,jk)
!                       buffwn (ind)= wnIO (ji,jj,jk)
!                       buffavt(ind)= avtIO(ji,jj,jk)
!                  enddo
!                 enddo
!                enddo
!
!                do jj =1 , jpj
!                 do ji =1 , jpi
!                   ind           = ji + jpi * (jj-1)
!                   buffvatm (ind)= vatmIO(ji,jj)
!                   buffemp  (ind)= empIO (ji,jj)
!                   buffqsr  (ind)= qsrIO (ji,jj)
!                   buffbblx (ind)= bblxIO(ji,jj)
!                   buffbbly (ind)= bblyIO(ji,jj)
!                  enddo
!                 enddo
!         ENDIF  ! if jn = 1

! ******  2.send buffer to rank 0
          call MPI_SEND(jpi, 1,mpi_integer, 0, 1, mpi_comm_world,ierr)
          call MPI_SEND(jpj, 1,mpi_integer, 0, 2, mpi_comm_world, ierr)
          call MPI_SEND(nimpp, 1,mpi_integer, 0, 3, mpi_comm_world,ierr)
          call MPI_SEND(njmpp, 1,mpi_integer, 0, 4, mpi_comm_world,ierr)
          call MPI_SEND(nlei, 1,mpi_integer, 0, 5, mpi_comm_world,ierr)
          call MPI_SEND(nlej, 1,mpi_integer, 0, 6, mpi_comm_world,ierr)
          call MPI_SEND(nldi, 1,mpi_integer, 0, 7, mpi_comm_world,ierr)
          call MPI_SEND(nldj, 1,mpi_integer, 0, 8, mpi_comm_world,ierr)
          call MPI_SEND(bufftrn, jpi * jpj * jpk,mpi_real8  , 0, 9, mpi_comm_world,ierr)
          call MPI_SEND(bufftrn_FN, jpi * jpj * jpk,mpi_real8  , 0, 10, mpi_comm_world,ierr)
!          IF (jn .EQ. 1) THEN
!               call MPI_SEND(buffsn  , jpi*jpj*jpk  ,mpi_real8, 0, 11, mpi_comm_world,ierr)
!               call MPI_SEND(bufftn  , jpi*jpj*jpk  ,mpi_real8, 0, 12, mpi_comm_world,ierr)
!               call MPI_SEND(buffvatm, jpi*jpj      ,mpi_real8, 0, 13, mpi_comm_world,ierr)
!               call MPI_SEND(buffemp , jpi*jpj      ,mpi_real8, 0, 14, mpi_comm_world,ierr)
!               call MPI_SEND(buffqsr , jpi*jpj      ,mpi_real8, 0, 15, mpi_comm_world,ierr)
!               call MPI_SEND(buffun  , jpi*jpj*jpk  ,mpi_real8, 0, 16, mpi_comm_world,ierr)
!               call MPI_SEND(buffbblx, jpi*jpj      ,mpi_real8, 0, 17, mpi_comm_world,ierr)
!               call MPI_SEND(buffvn  , jpi*jpj*jpk  ,mpi_real8, 0, 18, mpi_comm_world,ierr)
!               call MPI_SEND(buffbbly, jpi*jpj      ,mpi_real8, 0, 19, mpi_comm_world,ierr)
!               call MPI_SEND(buffwn  , jpi*jpj*jpk  ,mpi_real8, 0, 20, mpi_comm_world,ierr)
!               call MPI_SEND(buffavt , jpi*jpj*jpk  ,mpi_real8, 0, 21, mpi_comm_world,ierr)
!               call MPI_SEND(buffpp  , jpi*jpj*jpk*3,mpi_real8, 0, 22, mpi_comm_world,ierr)
!          ENDIF


      endif ! IF LABEL 1, if(rank == 0)

!************* END COLLECTING DATA  *****************



! *********** START WRITING **************************
      if(rank == 0) then ! IF LABEL 4

            DO ji=1, jpiglo ; lon(ji)= REAL(glamt(1,1),4)+(ji-1)*0.125 ; ENDDO
            DO jj=1, jpjglo ; lat(jj)= REAL(gphit(1,1),4)+(jj-1)*0.125 ; ENDDO
            DO jk=1, jpk    ; dep(jk)= REAL(gdept(jk),4)               ; ENDDO


            d2f3d= REAL(tottrnIO,4) !  conversion from double to float

         CALL IOOGSNC_AVE_3D2(output_file_nc, var_nc, jpiglo, jpjglo,jpk,lon,lat,dep, datefrom, dateTo, d2f3d)

            write(4445,rec=1) (((REAL(tottrnFN(ji,jj,jk),4), ji=1,jpiglo), jj=1, jpjglo), jk=1, jpk)

!            IF (jn .EQ. 1) THEN
!                write(4446,rec=1) (((REAL(totsnIO (ji,jj,jk),4), ji=1,jpiglo), jj=1,jpjglo), jk=1,jpk)
!                write(4446,rec=2) (((REAL(tottnIO (ji,jj,jk),4), ji=1,jpiglo), jj=1,jpjglo), jk=1,jpk)
!                write(4446,rec=3) (((REAL(totunIO (ji,jj,jk),4), ji=1,jpiglo), jj=1,jpjglo), jk=1,jpk)
!                write(4446,rec=4) (((REAL(totvnIO (ji,jj,jk),4), ji=1,jpiglo), jj=1,jpjglo), jk=1,jpk)
!                write(4446,rec=5) (((REAL(totwnIO (ji,jj,jk),4), ji=1,jpiglo), jj=1,jpjglo), jk=1,jpk)
!                write(4446,rec=6) (((REAL(totavtIO(ji,jj,jk),4), ji=1,jpiglo), jj=1,jpjglo), jk=1,jpk)
!                write(4447,rec=1) ((REAL(totvatmIO(ji,jj)   ,4), ji=1,jpiglo), jj=1,jpjglo)
!                write(4447,rec=2) ((REAL(totempIO (ji,jj)   ,4), ji=1,jpiglo), jj=1,jpjglo)
!                write(4447,rec=3) ((REAL(totqsrIO (ji,jj)   ,4), ji=1,jpiglo), jj=1,jpjglo)
!                write(4447,rec=4) ((REAL(totbblxIO(ji,jj)   ,4), ji=1,jpiglo), jj=1,jpjglo)
!                write(4447,rec=5) ((REAL(totbblyIO(ji,jj)   ,4), ji=1,jpiglo), jj=1,jpjglo)
!
!
! ******************  output GROSS  PRIMARY  PRODUCTION *******************
!
!               d2f3d= REAL(totppgIO,4)     !   conversion from double to float
!               var(1:3) = 'ppg'
!               CALL IOOGSNC_AVE_3D2(ppg_file_nc, var,jpiglo,jpjglo,jpk,lon,lat,dep,datefrom, dateTo,d2f3d)
!
! *****************  output NET  PRIMARY  PRODUCTION ********************
!
!               d2f3d= REAL(totppnIO,4)    !   conversion from double to float
!               var(1:3) = 'ppn'
!               CALL IOOGSNC_AVE_3D2(ppg_file_nc, var,jpiglo,jpjglo,jpk,lon,lat,dep,datefrom, dateTo,d2f3d)
!
********************* output Bacterial  PRODUCTION *********************
!               d2f3d= REAL(totppbIO,4)    !   conversion from double to float
!               var(1:3) = 'ppb'
!               CALL IOOGSNC_AVE_3D2(ppg_file_nc, var,jpiglo,jpjglo,jpk,lon,lat,dep,datefrom, dateTo,d2f3d)
!
!            ENDIF ! IF (jn .EQ. 1) THEN

            close(4445)
!            close(4446)
!            close(4447)


      end if ! IF LABEL 4  if(rank == 0)




      traIO(:,:,:,jn) = 0.  !      we reset matrix for new average
!      IF(jn .EQ. 1) THEN
!          snIO     = 0.
!          tnIO     = 0.
!          vatmIO   = 0.
!          empIO    = 0.
!          qsrIO    = 0.
!          unIO     = 0.
!          bblxIO   = 0.
!          vnIO     = 0.
!          bblyIO   = 0.
!          wnIO     = 0.
!          avtIO    = 0.
!          tra_ppIO = 0.
!      ENDIF

      END DO ! do jn=1,jptra , DO LABEL 5


      if(rank == 0) WRITE(numout,*) '**** trcdit : write NetCDF passive tracer concentration'


      END SUBROUTINE trcdit
