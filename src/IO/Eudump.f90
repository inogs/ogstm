      SUBROUTINE Eudump(datemean,datefrom,dateTo,FREQ_GROUP)
!     ******************
      USE calendar
      USE myalloc
      USE IO_mem
      USE FN_mem
      USE OPT_mem
      USE TIME_MANAGER
      use mpi
      USE ogstm_mpi_module

      IMPLICIT NONE


      CHARACTER(LEN=17), INTENT(IN) :: datemean, dateFrom, dateTo
      INTEGER, INTENT(IN) :: FREQ_GROUP
      INTEGER jk,jj,ji, jn, jn_high
      INTEGER ind, i_contribution, j_contribution
      CHARACTER(10) newfile
      CHARACTER(LEN=60) bkpname
      CHARACTER(LEN=11) DIR
      logical IsBackup
      integer ave_counter


      CHARACTER(LEN=56) dia_file_nc
      CHARACTER(LEN=20)  var
      character(len=4) :: wave_length 

      INTEGER idrank, ierr, istart, jstart, iPe, iPd, jPe, jPd, status(MPI_STATUS_SIZE)
      INTEGER irange, jrange
      INTEGER totistart, totiend, relistart, reliend
      INTEGER totjstart, totjend, reljstart, reljend
      double precision ::  Miss_val =1.e20



     ! call mppsync()
! ----------------------------------------
      IsBackup =  (datemean.eq.dateTo)
      if (lwp) write(*,*) 'diadump IsBackup = ',IsBackup, ' group ' ,FREQ_GROUP
! ----------------------------------------

      SELECT CASE (FREQ_GROUP)
        CASE (1) 
       ave_counter=ave_counter_1 
       DIR='AVE_FREQ_1/'
        CASE (2) 
       ave_counter=ave_counter_2 
       DIR='AVE_FREQ_2/'
      END SELECT

       jn_high = 0
! ! ******************  3D RADIATIVE TRANSPORT OUTPUT   *******************
       DO jn =1 , nlt

           if (.not.is_time_to_save(jn,FREQ_GROUP,3)) CYCLE
           if (FREQ_GROUP.eq.1) jn_high = jn_high+1

           if (myrank == 0) then                    ! IF LABEL 1


! ! ******* myrank 0 sets indexes of tot matrix where to place its own part

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

      if (FREQ_GROUP.eq.1) then
      tottrnIO (:,totjstart:totjend,totistart:totiend) = Eu_DIA_IO_HIGH(:,reljstart:reljend,relistart:reliend,jn_high)
      else
      tottrnIO (:, totjstart:totjend,totistart:totiend) = Eu_DIA_IO(:, reljstart:reljend, relistart:reliend,jn) ! diagnostic from reaction model
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
                 call MPI_RECV(buffDIA  ,jpi_rec*jpj_rec*jpk  ,mpi_real8   ,idrank, 9,mpi_comm_world, status, ierr)
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
                 do ji =totistart,totiend ! 3d vars
                    i_contribution = jpk*jpj_rec*(ji-1 -totistart+ relistart )
                    do jj =totjstart,totjend
                       j_contribution = jpk*(jj-totjstart+ reljstart-1)
                       do jk =1 , jpk
                            ind = jk + j_contribution + i_contribution
                            tottrnIO(jk,jj,ji)= buffDIA (ind)
                       enddo
                    enddo
                 enddo
              enddo !idrank = 1, size-1
           else  ! IF LABEL 1,  if(myrank == 0)
       if (FREQ_GROUP.eq.2) then
             do ji =1, jpi
               i_contribution= jpk*jpj * (ji - 1 )
               do jj =1 , jpj
                 j_contribution=jpk*(jj-1)
                 do jk =1 , jpk
                   ind         =  jk + j_contribution + i_contribution
                   buffDIA(ind) = Eu_DIA_IO(jk,jj,ji,jn)
                 enddo
               enddo
             enddo
       else
             do ji =1, jpi
               i_contribution= jpk*jpj * (ji - 1 )
               do jj =1 , jpj
                 j_contribution=jpk*(jj-1)
                 do jk =1 , jpk
                   ind         =  jk + j_contribution + i_contribution
                   buffDIA(ind) = Eu_DIA_IO_HIGH(jk,jj,ji,jn_high)
                 enddo
               enddo
             enddo
       endif
               call MPI_SEND(jpi  , 1,mpi_integer, 0, 1, mpi_comm_world,ierr)
               call MPI_SEND(jpj  , 1,mpi_integer, 0, 2, mpi_comm_world,ierr)
               call MPI_SEND(nimpp, 1,mpi_integer, 0, 3, mpi_comm_world,ierr)
               call MPI_SEND(njmpp, 1,mpi_integer, 0, 4, mpi_comm_world,ierr)
               call MPI_SEND(nlei , 1,mpi_integer, 0, 5, mpi_comm_world,ierr)
               call MPI_SEND(nlej , 1,mpi_integer, 0, 6, mpi_comm_world,ierr)
               call MPI_SEND(nldi , 1,mpi_integer, 0, 7, mpi_comm_world,ierr)
               call MPI_SEND(nldj , 1,mpi_integer, 0, 8, mpi_comm_world,ierr)
               call MPI_SEND(buffDIA  , jpk*jpj*jpi,mpi_real8, 0, 9, mpi_comm_world,ierr)
       endif ! IF LABEL 1, if(myrank == 0)
!************* END COLLECTING DATA  *****************

! *********** START WRITING **************************
      if (myrank == 0) then
           write(wave_length,"(I0.4)"), lam(jn)         
           var        =  'Eu_'//wave_length

          bkpname     = DIR//'ave.'//datemean//'.'//trim(var)//'.nc.bkp'
          dia_file_nc = DIR//'ave.'//datemean//'.'//trim(var)//'.nc'
              
          if (IsBackup) then
               CALL WRITE_AVE_BKP(bkpname,var,datefrom, dateTo,tottrnIO,ave_counter,deflate_ave, deflate_level_ave)
          else
               CALL WRITE_AVE(dia_file_nc,var,datefrom,dateTo, tottrnIO,deflate_ave, deflate_level_ave)
          endif


      end if ! IF LABEL 4  if(myrank == 0)
         if (.not.IsBackup) then
             if (FREQ_GROUP.eq.2) then
                Eu_DIA_IO(:,:,:,jn) = 0.
              else
                Eu_DIA_IO_HIGH(:,:,:,jn_high) = 0.
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
            IF (diahf(jn).eq.1)   IS_TIME_TO_SAVE = .true.
        ELSE
           IF (diahf_2d(jn).eq.1) IS_TIME_TO_SAVE = .true.
        ENDIF
      ENDIF
      END FUNCTION IS_TIME_TO_SAVE

      end SUBROUTINE Eudump
