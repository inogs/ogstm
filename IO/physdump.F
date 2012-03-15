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
      INTEGER ind

      CHARACTER*32 forcing2d_file, forcing3d_file


!      CHARACTER*29 ppg_file_nc, ppn_file_nc,ppb_file_nc! ave.20091231-12:00:00.P1n.nc
!      CHARACTER*29 dia1_file_nc, dia2_file_nc, dia3_file_nc, dia4_file_nc, dia5_file_nc
!      CHARACTER*29 dia6_file_nc, dia7_file_nc, dia8_file_nc, dia9_file_nc, dia10_file_nc
!      CHARACTER*29 dia11_file_nc, dia12_file_nc, dia13_file_nc, dia14_file_nc, dia15_file_nc
!      CHARACTER*29 dia16_file_nc, dia17_file_nc, dia18_file_nc, dia19_file_nc, dia20_file_nc
!      CHARACTER*29 dia21_file_nc, dia22_file_nc, dia23_file_nc, dia24_file_nc, dia25_file_nc
!      CHARACTER*29 dia26_file_nc, dia27_file_nc, dia28_file_nc, dia29_file_nc, dia30_file_nc
!      CHARACTER*29 dia31_file_nc, dia32_file_nc, dia33_file_nc

      CHARACTER*29 dia_file_nc
      CHARACTER*4  var, var_nc

      INTEGER :: ir, ir2d, ir3d
      INTEGER idrank, ierr, istart, jstart, iPe, iPd, jPe, jPd, status(MPI_STATUS_SIZE)
      INTEGER irange, jrange
      INTEGER totistart, totiend, relistart, reliend
      INTEGER totjstart, totjend, reljstart, reljend


      call mppsync()


!      ave_counter = 0   !  we reset the counter


      forcing3D_file = 'ave.'//datemean//'.grid3D.dat'
      forcing2D_file = 'ave.'//datemean//'.grid2D.dat'

!      dia1_file_nc    = 'ave.'//datemean//'.ppg.nc'//CHAR(0)
!      dia2_file_nc    = 'ave.'//datemean//'.ppn.nc'//CHAR(0)
!      dia3_file_nc    = 'ave.'//datemean//'.ppb.nc'//CHAR(0)
!      dia4_file_nc    = 'ave.'//datemean//'.DIC.nc'//CHAR(0)
!      dia5_file_nc    = 'ave.'//datemean//'.CO2.nc'//CHAR(0)
!      dia6_file_nc    = 'ave.'//datemean//'.HCO.nc'//CHAR(0)
!      dia7_file_nc    = 'ave.'//datemean//'.CO3.nc'//CHAR(0)
!      dia8_file_nc    = 'ave.'//datemean//'.AC_.nc'//CHAR(0)
!      dia9_file_nc    = 'ave.'//datemean//'.PH_.nc'//CHAR(0)
!      dia10_file_nc   = 'ave.'//datemean//'.pCO.nc'//CHAR(0)
!
!      dia11_file_nc   = 'ave.'//datemean//'.F01.nc'//CHAR(0)
!      dia12_file_nc   = 'ave.'//datemean//'.F02.nc'//CHAR(0)
!      dia13_file_nc   = 'ave.'//datemean//'.F03.nc'//CHAR(0)
!      dia14_file_nc   = 'ave.'//datemean//'.F04.nc'//CHAR(0)
!      dia15_file_nc   = 'ave.'//datemean//'.F05.nc'//CHAR(0)
!      dia16_file_nc   = 'ave.'//datemean//'.F06.nc'//CHAR(0)
!      dia17_file_nc   = 'ave.'//datemean//'.F07.nc'//CHAR(0)
!      dia18_file_nc   = 'ave.'//datemean//'.F08.nc'//CHAR(0)
!      dia19_file_nc   = 'ave.'//datemean//'.F09.nc'//CHAR(0)
!      dia20_file_nc   = 'ave.'//datemean//'.F10.nc'//CHAR(0)
!
!      dia21_file_nc   = 'ave.'//datemean//'.F11.nc'//CHAR(0)
!      dia22_file_nc   = 'ave.'//datemean//'.F12.nc'//CHAR(0)
!      dia23_file_nc   = 'ave.'//datemean//'.F13.nc'//CHAR(0)
!      dia24_file_nc   = 'ave.'//datemean//'.F14.nc'//CHAR(0)
!      dia25_file_nc   = 'ave.'//datemean//'.F15.nc'//CHAR(0)
!      dia26_file_nc   = 'ave.'//datemean//'.F16.nc'//CHAR(0)
!      dia27_file_nc   = 'ave.'//datemean//'.F17.nc'//CHAR(0)
!      dia28_file_nc   = 'ave.'//datemean//'.F18.nc'//CHAR(0)
!      dia29_file_nc   = 'ave.'//datemean//'.F19.nc'//CHAR(0)
!      dia30_file_nc   = 'ave.'//datemean//'.F20.nc'//CHAR(0)
!
!      dia31_file_nc   = 'ave.'//datemean//'.F21.nc'//CHAR(0)
!      dia32_file_nc   = 'ave.'//datemean//'.F22.nc'//CHAR(0)
!      dia33_file_nc   = 'ave.'//datemean//'.F23.nc'//CHAR(0)

      var   (4:4)    = CHAR(0)
      var_nc(4:4)    = CHAR(0)


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


! ******************  output GROSS  PRIMARY  PRODUCTION *******************
           do jn =1 , jptra_dia
              d2f3d= REAL(totdiaIO(:,:,:,jn),4)     !   conversion from double to float
              dia_file_nc ='ave.'//datemean//'.'//dia(jn)//'.nc'//CHAR(0)
              var(1:3) = dia(jn) ;
              CALL IOOGSNC_AVE_3D2(dia_file_nc, var,jpiglo,jpjglo,jpk,lon,lat,dep,datefrom, dateTo,d2f3d)
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
********************* output Bacterial  PRODUCTION *********************
!           d2f3d= REAL(totdiaIO(:,:,:,3),4)    !   conversion from double to float
!           var(1:3) = 'ppb'
!           CALL IOOGSNC_AVE_3D2(dia3_file_nc, var,jpiglo,jpjglo,jpk,lon,lat,dep,datefrom, dateTo,d2f3d)
! ******************  output DIC *******************
!
!           d2f3d= REAL(totdiaIO(:,:,:,4),4)     !   conversion from double to float
!           var(1:3) = 'DIC'
!           CALL IOOGSNC_AVE_3D2(dia4_file_nc, var,jpiglo,jpjglo,jpk,lon,lat,dep,datefrom, dateTo,d2f3d)
!
! *****************  output CO2 ********************
!
!           d2f3d= REAL(totdiaIO(:,:,:,5),4)    !   conversion from double to float
!           var(1:3) = 'CO2'
!           CALL IOOGSNC_AVE_3D2(dia5_file_nc, var,jpiglo,jpjglo,jpk,lon,lat,dep,datefrom, dateTo,d2f3d)
!
! ********************* output HCO3 *********************
!           d2f3d= REAL(totdiaIO(:,:,:,6),4)    !   conversion from double to float
!           var(1:3) = 'HCO'
!           CALL IOOGSNC_AVE_3D2(dia6_file_nc, var,jpiglo,jpjglo,jpk,lon,lat,dep,datefrom, dateTo,d2f3d)
!
! ******************  output CO3 *******************
!
!           d2f3d= REAL(totdiaIO(:,:,:,7),4)     !   conversion from double to float
!           var(1:3) = 'CO3'
!           CALL IOOGSNC_AVE_3D2(dia7_file_nc, var,jpiglo,jpjglo,jpk,lon,lat,dep,datefrom, dateTo,d2f3d)
!
! *****************  output AC ********************
!
!           d2f3d= REAL(totdiaIO(:,:,:,8),4)    !   conversion from double to float
!           var(1:3) = 'AC_'
!           CALL IOOGSNC_AVE_3D2(dia8_file_nc, var,jpiglo,jpjglo,jpk,lon,lat,dep,datefrom, dateTo,d2f3d)
!
********************* output Bacterial  HCO3 *********************
!           d2f3d= REAL(totdiaIO(:,:,:,9),4)    !   conversion from double to float
!           var(1:3) = 'PH_'
!           CALL IOOGSNC_AVE_3D2(dia9_file_nc, var,jpiglo,jpjglo,jpk,lon,lat,dep,datefrom, dateTo,d2f3d)
!
********************* pCO2 *********************
!
!           d2f3d= REAL(totdiaIO(:,:,:,10),4)    !   conversion from double to float
!           var(1:3) = 'pCO'
!           CALL IOOGSNC_AVE_3D2(dia10_file_nc, var,jpiglo,jpjglo,jpk,lon,lat,dep,datefrom, dateTo,d2f3d)
!
!
! ******************  output F01 *******************
!
!           d2f3d= REAL(totdiaIO(:,:,:,11),4)     !   conversion from double to float
!           var(1:3) = 'F01'
!           CALL IOOGSNC_AVE_3D2(dia11_file_nc, var,jpiglo,jpjglo,jpk,lon,lat,dep,datefrom, dateTo,d2f3d)
!
! *****************  ********************
!
!           d2f3d= REAL(totdiaIO(:,:,:,12),4)    !   conversion from double to float
!           var(1:3) = 'F02'
!           CALL IOOGSNC_AVE_3D2(dia12_file_nc, var,jpiglo,jpjglo,jpk,lon,lat,dep,datefrom, dateTo,d2f3d)
!
*********************  *********************
!           d2f3d= REAL(totdiaIO(:,:,:,13),4)    !   conversion from double to float
!           var(1:3) = 'F03'
!           CALL IOOGSNC_AVE_3D2(dia13_file_nc, var,jpiglo,jpjglo,jpk,lon,lat,dep,datefrom, dateTo,d2f3d)
! ******************   *******************
!
!           d2f3d= REAL(totdiaIO(:,:,:,14),4)     !   conversion from double to float
!           var(1:3) = 'F04'
!           CALL IOOGSNC_AVE_3D2(dia14_file_nc, var,jpiglo,jpjglo,jpk,lon,lat,dep,datefrom, dateTo,d2f3d)
!
! *****************   ********************
!
!           d2f3d= REAL(totdiaIO(:,:,:,15),4)    !   conversion from double to float
!           var(1:3) = 'F05'
!           CALL IOOGSNC_AVE_3D2(dia15_file_nc, var,jpiglo,jpjglo,jpk,lon,lat,dep,datefrom, dateTo,d2f3d)
!
*********************  *********************
!           d2f3d= REAL(totdiaIO(:,:,:,16),4)    !   conversion from double to float
!           var(1:3) = 'F06'
!           CALL IOOGSNC_AVE_3D2(dia16_file_nc, var,jpiglo,jpjglo,jpk,lon,lat,dep,datefrom, dateTo,d2f3d)
!
! ******************   *******************
!
!           d2f3d= REAL(totdiaIO(:,:,:,17),4)     !   conversion from double to float
!           var(1:3) = 'F07'
!           CALL IOOGSNC_AVE_3D2(dia17_file_nc, var,jpiglo,jpjglo,jpk,lon,lat,dep,datefrom, dateTo,d2f3d)
!
! *****************   ********************
!
!           d2f3d= REAL(totdiaIO(:,:,:,18),4)    !   conversion from double to float
!           var(1:3) = 'F08'
!           CALL IOOGSNC_AVE_3D2(dia18_file_nc, var,jpiglo,jpjglo,jpk,lon,lat,dep,datefrom, dateTo,d2f3d)
!
*********************  *********************
!           d2f3d= REAL(totdiaIO(:,:,:,19),4)    !   conversion from double to float
!           var(1:3) = 'F09'
!           CALL IOOGSNC_AVE_3D2(dia19_file_nc, var,jpiglo,jpjglo,jpk,lon,lat,dep,datefrom, dateTo,d2f3d)
!
*********************  *********************
!           d2f3d= REAL(totdiaIO(:,:,:,20),4)    !   conversion from double to float
!           var(1:3) = 'F10'
!           CALL IOOGSNC_AVE_3D2(dia20_file_nc, var,jpiglo,jpjglo,jpk,lon,lat,dep,datefrom, dateTo,d2f3d)
!
! ******************  output F21 *******************
!
!           d2f3d= REAL(totdiaIO(:,:,:,21),4)     !   conversion from double to float
!           var(1:3) = 'F11'
!           CALL IOOGSNC_AVE_3D2(dia21_file_nc, var,jpiglo,jpjglo,jpk,lon,lat,dep,datefrom, dateTo,d2f3d)
!
! *****************  ********************
!
!           d2f3d= REAL(totdiaIO(:,:,:,22),4)    !   conversion from double to float
!           var(1:3) = 'F12'
!           CALL IOOGSNC_AVE_3D2(dia22_file_nc, var,jpiglo,jpjglo,jpk,lon,lat,dep,datefrom, dateTo,d2f3d)
!
*********************  *********************
!           d2f3d= REAL(totdiaIO(:,:,:,23),4)    !   conversion from double to float
!           var(1:3) = 'F13'
!           CALL IOOGSNC_AVE_3D2(dia23_file_nc, var,jpiglo,jpjglo,jpk,lon,lat,dep,datefrom, dateTo,d2f3d)
! ******************   *******************
!
!           d2f3d= REAL(totdiaIO(:,:,:,24),4)     !   conversion from double to float
!           var(1:3) = 'F14'
!           CALL IOOGSNC_AVE_3D2(dia24_file_nc, var,jpiglo,jpjglo,jpk,lon,lat,dep,datefrom, dateTo,d2f3d)
!
 ! **   ********************
!
!           d2f3d= REAL(totdiaIO(:,:,:,25),4)    !   conversion from double to float
!           var(1:3) = 'F15'
!           CALL IOOGSNC_AVE_3D2(dia25_file_nc, var,jpiglo,jpjglo,jpk,lon,lat,dep,datefrom, dateTo,d2f3d)
!
! *********************  *********************
!           d2f3d= REAL(totdiaIO(:,:,:,26),4)    !   conversion from double to float
!           var(1:3) = 'F16'
!           CALL IOOGSNC_AVE_3D2(dia26_file_nc, var,jpiglo,jpjglo,jpk,lon,lat,dep,datefrom, dateTo,d2f3d)
!
 ! ***   *******************
!
!           d2f3d= REAL(totdiaIO(:,:,:,27),4)     !   conversion from double to float
!           var(1:3) = 'F17'
!           CALL IOOGSNC_AVE_3D2(dia27_file_nc, var,jpiglo,jpjglo,jpk,lon,lat,dep,datefrom, dateTo,d2f3d)
!
 ! **   ********************
!
!           d2f3d= REAL(totdiaIO(:,:,:,28),4)    !   conversion from double to float
!           var(1:3) = 'F18'
!           CALL IOOGSNC_AVE_3D2(dia28_file_nc, var,jpiglo,jpjglo,jpk,lon,lat,dep,datefrom, dateTo,d2f3d)
!
*********************  *********************
!           d2f3d= REAL(totdiaIO(:,:,:,29),4)    !   conversion from double to float
!           var(1:3) = 'F19'
!           CALL IOOGSNC_AVE_3D2(dia29_file_nc, var,jpiglo,jpjglo,jpk,lon,lat,dep,datefrom, dateTo,d2f3d)
!
*********************  *********************
!           d2f3d= REAL(totdiaIO(:,:,:,30),4)    !   conversion from double to float
!           var(1:3) = 'F20'
!           CALL IOOGSNC_AVE_3D2(dia30_file_nc, var,jpiglo,jpjglo,jpk,lon,lat,dep,datefrom, dateTo,d2f3d)
!
 ! ***  output F21 *******************
!
!           d2f3d= REAL(totdiaIO(:,:,:,31),4)     !   conversion from double to float
!           var(1:3) = 'F21'
!           CALL IOOGSNC_AVE_3D2(dia31_file_nc, var,jpiglo,jpjglo,jpk,lon,lat,dep,datefrom, dateTo,d2f3d)
!
 ! **  ********************
!
!           d2f3d= REAL(totdiaIO(:,:,:,32),4)    !   conversion from double to float
!           var(1:3) = 'F22'
!           CALL IOOGSNC_AVE_3D2(dia32_file_nc, var,jpiglo,jpjglo,jpk,lon,lat,dep,datefrom, dateTo,d2f3d)
!
*********************  *********************
!           d2f3d= REAL(totdiaIO(:,:,:,33),4)    !   conversion from double to float
!           var(1:3) = 'F23'
!           CALL IOOGSNC_AVE_3D2(dia33_file_nc, var,jpiglo,jpjglo,jpk,lon,lat,dep,datefrom, dateTo,d2f3d)


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
