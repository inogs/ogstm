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
       ! epascolo USE myalloc_mpp
       USE TIME_MANAGER
       USE IO_MEM , ONLY : ave_counter_1, ave_counter_2, existFilebkp

       IMPLICIT NONE

!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER jn, jn_high
      CHARACTER(LEN=37) filename
      CHARACTER(LEN=43) bkpname
      logical existFile
      logical bkp1hasbeenread,bkp2hasbeenread


! 0. initialisations
       bkpname  = 'AVE_FREQ_1/ave.20111231-15:30:00.N1p.nc.bkp'
       filename = 'RESTARTS/RST.20111231-15:30:00.N1p.nc'
       bkp1hasbeenread = .false.
       bkp2hasbeenread = .false.


      IF(lwp) THEN
          WRITE(numout,*) ' '
          WRITE(numout,*) ' *** trcrst beginning of restart for'
          WRITE(numout,*) ' passive tracer'
          WRITE(numout,*) '   with the time TimeStepStart : ',TimeStepStart
          WRITE(numout,*) ' '
      ENDIF

      jn_high=0
      DO jn=1, jptra  ! global loop on tracers to read restart


         filename = 'RESTARTS/RST.'//DateStart//'.'//trim(ctrcnm(jn))//'.nc'
         CALL readnc_slice_double(filename, 'TRN'//trim(ctrcnm(jn)), trn(:,:,:,jn) )



! ********************   we put initial undef to 0
          trn(:,:,:,jn) = trn(:,:,:,jn) * tmask
      
          trb = trn

!         SECTION AVE backup
!         rsttrn var is used instead of to define another one

          bkpname= 'AVE_FREQ_2/ave.'//DateStart//'.'//trim(ctrcnm(jn))//'.nc.bkp'

          INQUIRE(FILE=bkpname, EXIST=existFilebkp)

          if (existFilebkp) then
             if (lwp) write(*,*) 'reading ', bkpname
             CALL readnc_slice_double(bkpname,ctrcnm(jn), traIO(:,:,:,jn) )
             if (.not.bkp2hasbeenread) then
               call get_att_int( bkpname,'ave_counter', ave_counter_2)
               call get_att_char(bkpname,'DateStart'  , BKPdatefrom_2)
               bkp2hasbeenread=.true.
             endif
          else
             traIO(:,:,:,jn) = trn(:,:,:,jn)
          endif



         IF (ctr_hf(jn).eq.1)  THEN
           jn_high = jn_high + 1
           bkpname= 'AVE_FREQ_1/ave.'//DateStart//'.'//trim(ctrcnm(jn))//'.nc.bkp'
           INQUIRE(FILE=bkpname, EXIST=existFilebkp)

           if (existFilebkp) then
             if (lwp) write(*,*) 'reading ', bkpname
             CALL readnc_slice_double(bkpname,ctrcnm(jn), traIO_HIGH(:,:,:,jn_high) )
                if (.not.bkp1hasbeenread) then

                  call get_att_int( bkpname,'ave_counter', ave_counter_1)
                  call get_att_char(bkpname,'DateStart'  , BKPdatefrom_1)
                  bkp1hasbeenread=.true.
                endif
           else
              traIO_HIGH(:,:,:,jn_high) = 0.0
           endif
         ENDIF


      END DO


! ******************** 3D DIAGNOSTICS  ***********************************
      jn_high=0
      DO jn=1, jptra_dia

          bkpname    = 'AVE_FREQ_2/ave.'//DateStart//'.'//trim(dianm(jn))//'.nc.bkp'
          INQUIRE(FILE=bkpname, EXIST=existFile)

          if (existFile) then
             if (lwp) write(*,*) 'reading ', bkpname
             CALL readnc_slice_double(bkpname,dianm(jn), tra_DIA_IO(jn,:,:,:) )
             if (.not.bkp2hasbeenread) then
                call get_att_int( bkpname,'ave_counter', ave_counter_2)
                call get_att_char(bkpname,'DateStart'  , BKPdatefrom_2)
                bkp2hasbeenread=.true.
             endif
          else
             tra_DIA_IO(jn,:,:,:) = 0.0
          endif

          IF (diahf(jn).eq.1)  THEN
           jn_high = jn_high + 1
           bkpname= 'AVE_FREQ_1/ave.'//DateStart//'.'//trim(dianm(jn))//'.nc.bkp'
           INQUIRE(FILE=bkpname, EXIST=existFile)
           if (existFile) then
             if (lwp) write(*,*) 'reading ', bkpname
             CALL readnc_slice_double(bkpname,trim(dianm(jn)), tra_DIA_IO_HIGH(:,:,:,jn_high) )
                    if (.not.bkp1hasbeenread) then
                      call get_att_int( bkpname,'ave_counter', ave_counter_1)
                      call get_att_char(bkpname,'DateStart'  , BKPdatefrom_1)
                      bkp1hasbeenread=.true.
                    endif
           else
              tra_DIA_IO_HIGH(:,:,:,jn_high) = 0.0
           endif

          ENDIF
      END DO


! ******************** 2D DIAGNOSTICS  ***********************************

      jn_high=0
      DO jn=1, jptra_dia_2d

          bkpname    = 'AVE_FREQ_2/ave.'//DateStart//'.'//trim(dianm_2d(jn))//'.nc.bkp'
          INQUIRE(FILE=bkpname, EXIST=existFile)

          if (existFile) then
             if (lwp) write(*,*) 'reading ', bkpname
             CALL readnc_slice_double(bkpname,dianm_2d(jn), tra_DIA_2d_IO(:,:,jn) )
             if (.not.bkp2hasbeenread) then
                call get_att_int( bkpname,'ave_counter', ave_counter_2)
                call get_att_char(bkpname,'DateStart'  , BKPdatefrom_2)
                bkp2hasbeenread=.true.
             endif
          else
             tra_DIA_2d_IO(:,:,jn) = 0.0
          endif

          IF (diahf_2d(jn).eq.1)  THEN
           jn_high = jn_high + 1
           bkpname= 'AVE_FREQ_1/ave.'//DateStart//'.'//trim(dianm_2d(jn))//'.nc.bkp'
           INQUIRE(FILE=bkpname, EXIST=existFile)
           if (existFile) then
             if (lwp) write(*,*) 'reading ', bkpname
             CALL readnc_slice_double(bkpname,trim(dianm(jn)), tra_DIA_2d_IO_HIGH(:,:,jn_high) )
                    if (.not.bkp1hasbeenread) then
                      call get_att_int( bkpname,'ave_counter', ave_counter_1)
                      call get_att_char(bkpname,'DateStart'  , BKPdatefrom_1)
                      bkp1hasbeenread=.true.
                    endif
           else
              tra_DIA_2d_IO_HIGH(:,:,jn_high) = 0.0
           endif

          ENDIF



      END DO





      bkpname    = 'AVE_PHYS/ave.'//DateStart//'.phys.nc.bkp'
      INQUIRE(FILE=bkpname, EXIST=existFilebkp)
      if (existFilebkp) then
          call readnc_slice_double(   bkpname,'vosaline',  snIO)
          call readnc_slice_double(   bkpname,'votemper',  tnIO)

          call readnc_slice_double(   bkpname,'vozocrtx',  unIO)
          call readnc_slice_double(   bkpname,'vomecrty',  vnIO)
          call readnc_slice_double(   bkpname,'vovecrtz',  wnIO)
          call readnc_slice_double(   bkpname,'votkeavt', avtIO)
          call readnc_slice_double(   bkpname,'e3t', e3tIO)

          call readnc_slice_double_2d(bkpname,'soshfldo', qsrIO)
          call readnc_slice_double_2d(bkpname,'sowindsp',vatmIO)
          call readnc_slice_double_2d(bkpname,'sowaflcd', empIO)
      else
          snIO      = 0.0
          tnIO      = 0.0
          vatmIO    = 0.0
          empIO     = 0.0
          qsrIO     = 0.0
          unIO      = 0.0
          vnIO      = 0.0
          wnIO      = 0.0
          avtIO     = 0.0
          e3tIO     = 0.0
      endif

      if (lwp) write(*,*) 'trcrst avecounter_1 = ', ave_counter_1, ' avecounter_2 = ', ave_counter_2
      IsStartBackup_1 = bkp1hasbeenread
      IsStartBackup_2 = bkp2hasbeenread


      END SUBROUTINE trcrst

