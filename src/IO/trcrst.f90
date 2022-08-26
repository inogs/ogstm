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
       USE TIME_MANAGER
       USE IO_MEM , ONLY : elapsed_time_1, elapsed_time_2, existFilebkp

       IMPLICIT NONE

!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER jn, jn_high
      CHARACTER(LEN=37) filename
      CHARACTER(LEN=100) bkpname
      logical existFile
      logical bkp1hasbeenread,bkp2hasbeenread
      CHARACTER(LEN=10) dirct

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


         filename = 'RESTARTS/RST.'//DateStart//'.'//trim(ctrcnm(jn))// & 
                '.nc'
         CALL read_double(filename, 'TRN'//trim(ctrcnm(jn)),0, trn(:,:,:,jn) )



! ********************   we put initial undef to 0
          trn(:,:,:,jn) = trn(:,:,:,jn) * tmask
      
          trb = trn

!         SECTION AVE backup
!         rsttrn var is used instead of to define another one

          bkpname= 'AVE_FREQ_2/ave.'//DateStart//'.'//trim(ctrcnm(jn))//'.nc.bkp'

          INQUIRE(FILE=bkpname, EXIST=existFilebkp)

          if (existFilebkp) then
             if (lwp) write(*,*) 'reading ', bkpname
             CALL read_double(bkpname,ctrcnm(jn),0, traIO(:,:,:,jn) )
             if (.not.bkp2hasbeenread) then
               call readnc_scalar_double(bkpname,'elapsed_time',elapsed_time_2)
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
             CALL read_double(bkpname,ctrcnm(jn),0, traIO_HIGH(:,:,:,jn_high) )
                if (.not.bkp1hasbeenread) then
                  call readnc_scalar_double(bkpname,'elapsed_time',elapsed_time_1)
                  call get_att_char(bkpname,'DateStart'  , BKPdatefrom_1)
                  bkp1hasbeenread=.true.
                endif
           else
              traIO_HIGH(:,:,:,jn_high) = trn(:,:,:,jn)
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
             CALL read_double(bkpname,dianm(jn),0, tra_DIA_IO(jn,:,:,:) )
             if (.not.bkp2hasbeenread) then
                call readnc_scalar_double(bkpname,'elapsed_time',elapsed_time_2)
                call get_att_char(bkpname,'DateStart'  , BKPdatefrom_2)
                bkp2hasbeenread=.true.
             endif
          else
             tra_DIA_IO(jn,:,:,:) = 0.0
          endif

          IF ((diahf(jn).eq.1).and.(diaWR(jn).eq.1))  THEN
           jn_high = jn_high + 1
           bkpname= 'AVE_FREQ_1/ave.'//DateStart//'.'//trim(dianm(jn))//'.nc.bkp'
           INQUIRE(FILE=bkpname, EXIST=existFile)
           if (existFile) then
             if (lwp) write(*,*) 'reading ', bkpname
             CALL read_double(bkpname,trim(dianm(jn)),0, tra_DIA_IO_HIGH(jn_high,:,:,:) )
                    if (.not.bkp1hasbeenread) then
                      call readnc_scalar_double(bkpname,'elapsed_time',elapsed_time_1)
                      call get_att_char(bkpname,'DateStart'  , BKPdatefrom_1)
                      bkp1hasbeenread=.true.
                    endif
           else
              tra_DIA_IO_HIGH(jn_high,:,:,:) = 0.0
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
             CALL read_double_2d(bkpname,dianm_2d(jn),0, tra_DIA_2d_IO(jn,:,:) )
             if (.not.bkp2hasbeenread) then
                call readnc_scalar_double(bkpname,'elapsed_time',elapsed_time_2)
                call get_att_char(bkpname,'DateStart'  , BKPdatefrom_2)
                bkp2hasbeenread=.true.
             endif
          else
             tra_DIA_2d_IO(jn,:,:) = 0.0
          endif

          IF ((diahf_2d(jn).eq.1).and.(diaWR_2d(jn).eq.1))  THEN
           jn_high = jn_high + 1
           bkpname= 'AVE_FREQ_1/ave.'//DateStart//'.'//trim(dianm_2d(jn))//'.nc.bkp'
           INQUIRE(FILE=bkpname, EXIST=existFile)
           if (existFile) then
             if (lwp) write(*,*) 'reading ', bkpname
             CALL read_double_2d(bkpname,trim(dianm_2d(jn)),0, tra_DIA_2d_IO_HIGH(jn_high,:,:) )
                    if (.not.bkp1hasbeenread) then
                      call readnc_scalar_double(bkpname,'elapsed_time',elapsed_time_1)
                      call get_att_char(bkpname,'DateStart'  , BKPdatefrom_1)
                      bkp1hasbeenread=.true.
                    endif
           else
              tra_DIA_2d_IO_HIGH(jn_high,:,:) = 0.0
           endif

          ENDIF



      END DO




      if(freq_ave_phys .eq. 1) then
        dirct="AVE_FREQ_1"
      else
        dirct="AVE_FREQ_2"
      end if
      
      bkpname = dirct//'/ave.'//DateStart//'.'//'vosaline'//'.nc.bkp'
        INQUIRE(FILE=bkpname, EXIST=existFilebkp)
        if (existFilebkp) then
          if(lwp)  write(*,*) 'reading ', bkpname
          call read_double(   bkpname,'vosaline',0,  snIO) 
        else 
          snIO = 0.0
        end if



      bkpname = dirct//'/ave.'//DateStart//'.'//'votemper'//'.nc.bkp'
        INQUIRE(FILE=bkpname, EXIST=existFilebkp)
        if (existFilebkp) then
          if(lwp)  write(*,*) 'reading ', bkpname
          call read_double(   bkpname,'votemper',0,  tnIO)      
        else
          tnIO = 0.0
        end if
     

      bkpname = dirct//'/ave.'//DateStart//'.'//'vozocrtx'//'.nc.bkp'
        INQUIRE(FILE=bkpname, EXIST=existFilebkp)
        if (existFilebkp) then
          if(lwp)  write(*,*) 'reading ', bkpname
          call read_double(   bkpname,'vozocrtx',0,  unIO)  
        else
          unIO = 0.0
        end if


      bkpname = dirct//'/ave.'//DateStart//'.'//'vomecrty'//'.nc.bkp'
        INQUIRE(FILE=bkpname, EXIST=existFilebkp)
        if (existFilebkp) then
          if(lwp)  write(*,*) 'reading ', bkpname
          call read_double(   bkpname,'vomecrty',0,  vnIO)  
        else
          vnIO = 0.0
        end if


      bkpname = dirct//'/ave.'//DateStart//'.'//'vovecrtz'//'.nc.bkp'
        INQUIRE(FILE=bkpname, EXIST=existFilebkp)
        if (existFilebkp) then
          if(lwp)  write(*,*) 'reading ', bkpname
          call read_double(   bkpname,'vovecrtz',0,  wnIO)  
        else
          wnIO = 0.0
        end if


      bkpname = dirct//'/ave.'//DateStart//'.'//'votkeavt'//'.nc.bkp'
        INQUIRE(FILE=bkpname, EXIST=existFilebkp)
        if (existFilebkp) then
          if(lwp)  write(*,*) 'reading ', bkpname
          call read_double(   bkpname,'votkeavt',0,  avtIO)  
        else
          avtIO = 0.0
        end if


      bkpname = dirct//'/ave.'//DateStart//'.'//'e3t'//'.nc.bkp'
        INQUIRE(FILE=bkpname, EXIST=existFilebkp)
        if (existFilebkp) then
          if(lwp)  write(*,*) 'reading ', bkpname
          call read_double(   bkpname,'e3t',0,  e3tIO)  
        else
          e3tIO = 0.0
        end if


      bkpname = dirct//'/ave.'//DateStart//'.'//'soshfldo'//'.nc.bkp'
        INQUIRE(FILE=bkpname, EXIST=existFilebkp)
        if (existFilebkp) then
          if(lwp)  write(*,*) 'reading ', bkpname
          call read_double_2d(   bkpname,'soshfldo',0, qsrIO)  
        else
          qsrIO = 0.0
        end if


      bkpname = dirct//'/ave.'//DateStart//'.'//'sowindsp'//'.nc.bkp'
        INQUIRE(FILE=bkpname, EXIST=existFilebkp)
        if (existFilebkp) then
          if(lwp)  write(*,*) 'reading ', bkpname
          call read_double_2d(   bkpname,'sowindsp',0,  vatmIO)  
        else
          vatmIO = 0.0
        end if


      bkpname = dirct//'/ave.'//DateStart//'.'//'sowaflcd'//'.nc.bkp'
        INQUIRE(FILE=bkpname, EXIST=existFilebkp)
        if (existFilebkp) then
          if(lwp)  write(*,*) 'reading ', bkpname
          call read_double_2d(   bkpname,'sowaflcd',0,  empIO)  
        else
          empIO = 0.0
        end if




      if (lwp) write(*,*) 'trcrst elapsed_time_1 = ', elapsed_time_1, &
                ' elapsed_time_2 = ', elapsed_time_2
      IsStartBackup_1 = bkp1hasbeenread
      IsStartBackup_2 = bkp2hasbeenread


      END SUBROUTINE trcrst

