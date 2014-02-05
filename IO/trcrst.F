      SUBROUTINE trcrst
!---------------------------------------------------------------------
!
!                       ROUTINE trcrst
!                     ******************
!
!  PURPOSE :
!  ---------
!     READ files for restart for passive tracer
!
!----------------------------------------------------------------------


       USE calendar
       USE myalloc
       USE myalloc_mpp
       USE TIME_MANAGER
       USE IO_MEM , ONLY : ave_counter,existFilebkp

       IMPLICIT NONE

!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER jn
      CHARACTER(LEN=28) filename
      CHARACTER(LEN=6)  varname
      CHARACTER(LEN=33) bkpname
      logical existFile

! 0. initialisations
       bkpname  = 'ave.20111231-15:30:00.N1p.nc.bkp'
       filename = 'RST.20111231-15:30:00.N1p.nc'
       varname  = 'TRNN1p'

!      qui leggeva il date.step


      IF(lwp) THEN
          WRITE(numout,*) ' '
          WRITE(numout,*) ' *** trcrst beginning of restart for'
          WRITE(numout,*) ' passive tracer'
          WRITE(numout,*) '   with the time TimeStepStart : ',TimeStepStart
          WRITE(numout,*) ' '
      ENDIF

      DO jn=1, jptra  ! global loop on tracers to read restart


         varname  = 'TRN'//ctrcnm(jn)
         filename = 'RST.'//DateStart//'.'//ctrcnm(jn)//'.nc'

         CALL readnc_slice_double(filename,varname, trn(:,:,:,jn) )

         varname = 'TRB'//ctrcnm(jn)
         CALL readnc_slice_double(filename,varname, trb(:,:,:,jn) )



! ********************   we put initial undef to 0
          trb(:,:,:,jn) = trb(:,:,:,jn) * tmask;
          trn(:,:,:,jn) = trn(:,:,:,jn) * tmask;

!         SECTION AVE backup
!         rsttrn var is used instead of to define another one

          bkpname= 'ave.'//DateStart//'.'//ctrcnm(jn)//'.nc.bkp'

          INQUIRE(FILE=bkpname, EXIST=existFilebkp)

          if (existFilebkp) then
              isStartBackup = .true.
              CALL readnc_slice_double(bkpname,ctrcnm(jn), traIO(:,:,:,jn) )


             call get_att_int(bkpname,'ave_counter', ave_counter)
             call get_att_char(bkpname,'DateStart', BKPdatefrom)



          endif

      END DO

      DO jn=1, jptra_dia

          bkpname    = 'ave.'//DateStart//'.'//dia(jn)//'.nc.bkp'
          INQUIRE(FILE=bkpname, EXIST=existFile)

          if (existFile) then
               CALL readnc_slice_double(bkpname,dia(jn), tra_ppIO(:,:,:,jn) )
          endif

      END DO

      bkpname    = 'ave.'//DateStart//'.phys.nc.bkp'
      if (existFilebkp) then
          call readnc_slice_double(bkpname,'vosaline',snIO);
          call readnc_slice_double(bkpname,'votemper',tnIO);

          call readnc_slice_double(bkpname,'vozocrtx',unIO);
          call readnc_slice_double(bkpname,'vomecrty',vnIO);
          call readnc_slice_double(bkpname,'vovecrtz',wnIO);
          call readnc_slice_double(bkpname,'votkeavt',avtIO);

          call readnc_slice_double_2d(bkpname,'soshfldo',qsrIO)
          call readnc_slice_double_2d(bkpname,'sowindsp',vatmIO)
          call readnc_slice_double_2d(bkpname,'sowaflcd',empIO)
      endif

      END SUBROUTINE trcrst

