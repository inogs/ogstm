      SUBROUTINE trcditDA(datemean,datefrom,dateTo)
!---------------------------------------------------------------------
!
!                       ROUTINE trcdit
!                     ******************
!
!  Purpose :
!  ---------
!     Standard output of passive tracer : concentration fields


      USE netcdf
      USE calendar
      USE myalloc
      USE IO_mem
      USE FN_mem
      USE TIME_MANAGER
      USE DA_mem
#ifdef key_mpp
      USE myalloc_mpp
#endif

      IMPLICIT NONE

! local declarations
! ==================
      CHARACTER(LEN=17), INTENT(IN) :: datemean, dateFrom, dateTo



      INTEGER ji, jj, jk, jn, jn_high
      INTEGER ind,jn2
      integer s, nc, counter
      integer timid, depid, yid, xid, idvar

      CHARACTER(LEN=39) output_file_nc  ! DA/ave.20091231-12:00:00.P1n.nc
      CHARACTER(LEN=3)  var

!----------------------------------------------------------------------
! statement functions
! ===================


      INTEGER idrank, ierr, istart, jstart, iPe, iPd, jPe, jPd, status(MPI_STATUS_SIZE)
      INTEGER irange, jrange
      INTEGER totistart, totiend, relistart, reliend
      INTEGER totjstart, totjend, reljstart, reljend

! ----------------------------------------

      call mppsync()

      if (lwp)  CHLtot = 0.0


!  Ghost Shells - Manual s Indexes-

      DO jn=1,17 ! DO LABEL 5
          var = varlistDA(jn)

          jn2 = GET_INDEX(var)
          jn_high = GET_HIGHFREQ(jn2)

          output_file_nc = 'DA__FREQ_1/ave.'//datemean//'.'//var//'.nc'

!         In order to avoid zeros in output_file_nc, when we reduce to float32
          do jk=1,jpk
          do jj=1,jpj
          do ji=1,jpi
           if (traIO_HIGH(ji,jj,jk,jn_high).lt.SMALL)  traIO_HIGH(ji,jj,jk,jn_high)=SMALL
          enddo
          enddo
          enddo

! *************** START COLLECTING DATA *****************************
      if(myrank == 0) then                    ! IF LABEL 1
          write(*,*) var, ' jn =',jn,jn2,jn_high
! ******* myrank 0 sets indexes of tot matrix where to place its own part
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


!***** START ASSEMBLING ***  myrank 0 puts its tracer part in the tot matrix

          tottrnIO(totistart:totiend, totjstart:totjend,:)= traIO_HIGH(relistart:reliend, reljstart:reljend,:,jn_high)



          do idrank = 1, size-1
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
              do jk =1 , jpk
               do jj =totjstart,totjend
                 do ji =totistart,totiend
                     ind = (ji-totistart+ relistart)+ (jj-totjstart+ reljstart-1)*jpi_rec+(jk-1)*jpj_rec*jpi_rec
                     tottrnIO(ji,jj,jk)= bufftrn   (ind)
                 enddo
                enddo
               enddo

          enddo !idrank = 1, size-1


      else  ! IF LABEL 1,  if(myrank == 0)
! **** work of the other myranks
! ****** 1. load  inf buffer their IO matrices


           do jk =1 , jpk
            do jj =1 , jpj
             do ji =1 , jpi
               ind            =  ji + jpi * (jj-1) + jpi * jpj *(jk-1)
               bufftrn   (ind)= traIO_HIGH( ji,jj,jk,jn_high)
             enddo
            enddo
          enddo





! ******  2.send buffer to myrank 0
          call MPI_SEND(jpi  , 1,mpi_integer, 0, 1, mpi_comm_world,ierr)
          call MPI_SEND(jpj  , 1,mpi_integer, 0, 2, mpi_comm_world,ierr)
          call MPI_SEND(nimpp, 1,mpi_integer, 0, 3, mpi_comm_world,ierr)
          call MPI_SEND(njmpp, 1,mpi_integer, 0, 4, mpi_comm_world,ierr)
          call MPI_SEND(nlei , 1,mpi_integer, 0, 5, mpi_comm_world,ierr)
          call MPI_SEND(nlej , 1,mpi_integer, 0, 6, mpi_comm_world,ierr)
          call MPI_SEND(nldi , 1,mpi_integer, 0, 7, mpi_comm_world,ierr)
          call MPI_SEND(nldj , 1,mpi_integer, 0, 8, mpi_comm_world,ierr)
          call MPI_SEND(bufftrn, jpi * jpj * jpk,mpi_real8  , 0, 9, mpi_comm_world,ierr)



      endif ! IF LABEL 1, if(myrank == 0)

!************* END COLLECTING DATA  *****************



! *********** START WRITING **************************
      if(myrank == 0) then ! IF LABEL 4

          d2f3d = REAL(tottrnIO,4)
          !write(*,*) 'writing ', output_file_nc
          CALL WRITE_AVE(output_file_nc,datefrom, dateTo, d2f3d);
          if (isaCHLVAR(VAR)) CHLtot = CHLtot + d2f3d ! remains in RAM, used by snutell

      end if ! IF LABEL 4  if(myrank == 0)


!         NIENTE AZZERAMENTO MATRICI

      END DO ! do jn=1,17 , DO LABEL 5


      if (myrank==0) then

        s = nf90_create(CHLSUP_FOR_DA, NF90_CLOBBER, nc)
        s= nf90_def_dim(nc,'lon'           , jpiglo,  xid)
        s= nf90_def_dim(nc,'lat'           , jpjglo,  yid)
        s= nf90_def_dim(nc,'depth'         , 1   ,depid)
        s= nf90_def_dim(nc,'time'  , NF90_UNLIMITED,timid)
        s = nf90_def_var(nc, 'lchlm', nf90_float, (/xid,yid,depid,timid/),  idVAR)

        s = nf90_put_att(nc,idVAR, 'missing_value',4.e+20)

        s =nf90_enddef(nc)

        CHL_SUP = CHLtot(:,:,1)
        s = nf90_put_var(nc, idVAR  ,CHL_SUP ); call handle_err1(s,counter,CHLSUP_FOR_DA)
        s=nf90_close(nc)


      endif







      CONTAINS
      INTEGER FUNCTION GET_INDEX(VARNAME)
      CHARACTER(LEN=3), INTENT(IN) :: VARNAME
      INTEGER i
      GET_INDEX=0
      do i =1, jptra
         if (ctrcnm(i).eq.VARNAME) GET_INDEX = i
      enddo
      END FUNCTION GET_INDEX

      INTEGER FUNCTION GET_HIGHFREQ(varind)
      integer, INTENT(IN) :: varind
      integer i
      GET_HIGHFREQ = 0
      do i=1,jptra_high
         if (highfreq_table(i).eq.varind) GET_HIGHFREQ = i
      enddo
      END FUNCTION GET_HIGHFREQ


      END SUBROUTINE trcditDA





