      SUBROUTINE fluxdump(datemean, datefrom, dateend,FREQ_GROUP)

!     ******************
!     Works only if file Fluxes.nc exists and for FREQ_GROUP=1


      USE DIA_mem
      USE netcdf
      use mpi
      USE IO_mem, only: elapsed_time_1, elapsed_time_2

      IMPLICIT NONE


      CHARACTER(LEN=17), INTENT(IN) :: datemean, datefrom, dateend
      INTEGER, INTENT(IN) :: FREQ_GROUP ! 1 = HIGH FREQ, 2 = LOW FREQ
      INTEGER jf,js,jn,jn_high, ji,jj,counter,counterV
      INTEGER idrank, ierr, istart,status(MPI_STATUS_SIZE)

      CHARACTER(LEN=32) flux_file


      integer s, nc,INDid
      integer nid, tid
      INTEGER IDS(jptra)
      double precision  Realcounter

      if (.not.existFileFluxes) RETURN
      IF  (freq_flux_dump.eq.1) THEN
           if (FREQ_GROUP.eq.2 )     RETURN
           Realcounter   =    1./elapsed_time_1
      ELSE
           if (FREQ_GROUP.eq.1 )     RETURN
           Realcounter   =    1./elapsed_time_2
      ENDIF

      flx_partTime = MPI_WTIME()

!!!!!!!!!!!!!!!!!!PROCESSOR ZERO COLLECTS ALL FLUXES AND DUMPS  ON NETCDF FILE
      if (myrank ==0) then
!     Rank Zero only creates netcdf file
!     FLUXES/flux.20021216-12:00:00.nc
          flux_file = 'FLUXES/flux.'//datemean//'.nc'


          ! Just to try without 'or'
          ! s = nf90_create(flux_file, or(nf90_clobber,NF90_HDF5), nc)
          s = nf90_create(flux_file, NF90_HDF5, nc)

          s = nf90_put_att(nc, nf90_global, 'Time_Start'     , datefrom)
          s = nf90_put_att(nc, nf90_global, 'Time___End'     ,  dateend)


          s= nf90_def_dim(nc,'n'           , FsizeGlo,  nid)
          s= nf90_def_dim(nc,'Type'        ,        7,  tid)

          ! ******************   var definition
          s = nf90_def_var(nc,'index' ,   nf90_int, (/nid/),   IndID)

          DO jn=1,jptra
             s = nf90_def_var(nc,ctrcnm(jn) ,   nf90_double, (/tid,nid/),   IDS(jn))
          ENDDO

           s =nf90_enddef(nc)
           counter = 1
           INDflxBuff = 0
           if (Fsize .GT. 0) THEN
               do jf=1,Fsize
                   ji = flx_ridxt(jf,4)
                   jj = flx_ridxt(jf,3)
                   INDflxDUMPZERO(jf)=INDflxDUMP(jf)
                   if ( (ji .EQ. 1) .OR. (ji .EQ. jpi) ) INDflxDUMPZERO(jf) = 0 ! Ghost cell value has index 0
                   if ( (jj .EQ. 1) .OR. (jj .EQ. jpj) ) INDflxDUMPZERO(jf) = 0 ! Ghost cell value has index 0
                   if (INDflxDUMPZERO(jf) .NE. 0) then
                       INDflxDUMPglo(counter) = INDflxDUMPZERO(jf)
                       counter = counter +1
                   endif
               enddo
           endif
           do idrank = 1,mpi_glcomm_size-1
               call MPI_RECV(INDflxBuff    , FsizeMax,                 mpi_integer, idrank, 1,mpi_comm_world, status, ierr)
               do jf=1,FsizeMax
                   if (INDflxBuff(jf) .NE. 0) then
                       INDflxDUMPglo(counter) = INDflxBuff(jf)
                       counter = counter +1
                   endif
               enddo
           end do
 !     index part
           !************************ now, put var
           s = nf90_put_var(nc, IndID,      INDflxDUMPglo)

           DO jn=1,jptra

               MflxDumpGlo = 0
               counterV=1
               if (Fsize .GT. 0) THEN
                   do jf=1,Fsize
                       if (INDflxDUMPZERO(jf) .NE. 0) then
                           do js=1,7
                               MflxDumpGlo(js,counterV) = diaflx(js,jf,jn)
                           end do
                           counterV = counterV +1
                       endif
                   enddo
               endif
               DO idrank = 1,mpi_glcomm_size-1
                   call MPI_RECV(INDflxBuff    , FsizeMax,    mpi_integer, idrank, 2,mpi_comm_world, status, ierr)
                   call MPI_RECV(diaflxBuff    , FsizeMax*7,  mpi_real8,   idrank, 3,mpi_comm_world, status, ierr)
                   DO jf=1,FsizeMax
                       if (INDflxBuff(jf) .NE. 0) then
                           DO js=1,7
                               MflxDumpGlo(js,counterV) = diaflxBuff(jf,js)
                           END DO
                           counterV = counterV +1
                       endif
                   END DO
               ENDDO ! loop on myrank for each tracers

      !************************ now, put var

          MflxDumpGlo = MflxDumpGlo * Realcounter ! we store average flux
          s = nf90_put_var(nc,IDS(jn) ,  MflxDumpGlo)
!         MflxDumpGlo_FLOAT = real(MflxDumpGlo,4)
!         s = nf90_put_var(nc,IDS(jn) ,  MflxDumpGlo_FLOAT)

          ENDDO ! loop on tracers


          s= nf90_close(nc)

!!!!!!!!!!!!!!!!!!END OF TASK ZERO JOB
      else ! Other ranks than zero
!!!!!!!!!!!!!!!!!!OTHER THAN RANK ZERO JOB STARTS --> SENDING DATA TO PROCESSOR ZERO

           INDflxBuff = 0
           diaflxBuff = 0
           if (Fsize .GT. 0) THEN
               do jf=1,Fsize
                   ji = flx_ridxt(jf,4)
                   jj = flx_ridxt(jf,3)
                   INDflxBuff(jf)=INDflxDUMP(jf)
                   if ( (ji .EQ. 1) .OR. (ji .EQ. jpi) ) INDflxBuff(jf) = 0 ! Ghost cell value has index 0
                   if ( (jj .EQ. 1) .OR. (jj .EQ. jpj) ) INDflxBuff(jf) = 0 ! Ghost cell value has index 0
               enddo
           endif
           call MPI_SEND(INDflxBuff, FsizeMax, mpi_integer, 0, 1, mpi_comm_world, ierr)
           DO jn=1,jptra
               IF (Fsize .GT. 0) THEN
                   DO jf=1,Fsize
                       if (INDflxBuff(jf) .NE. 0) then
                           do js=1,7
                               diaflxBuff(jf,js) = diaflx(js,jf,jn)
                           end do
                       endif
                   ENDDO
               ENDIF
               call MPI_SEND(INDflxBuff, FsizeMax, mpi_integer, 0, 2, mpi_comm_world, ierr)
               call MPI_SEND(diaflxBuff, FsizeMax*7, mpi_real8, 0, 3, mpi_comm_world, ierr)
           END DO ! loop on myrank for each tracers


      endif

!!!!!!!!!!!!!!!!!!END OF COMMUNICATION PART AND FILE DUMP

! reset integral of flux

      IF (Fsize.NE.0 ) diaflx = 0


      flx_partTime = MPI_WTIME() - flx_partTime
      flx_TotTime  = flx_TotTime + flx_partTime

      end SUBROUTINE fluxdump
