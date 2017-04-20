      SUBROUTINE mainASSIMILATION(datestr,datefrom)

      use filenames
      use DA_mem
      use mpi_str, only: Var3DCommunicator

      IMPLICIT NONE

      CHARACTER(LEN=17), INTENT(IN) :: datestr, dateFrom

      character(LEN=1024) SATFILE, VARFILE
      character(LEN=46) SUFFIX
      character(LEN=2) MONTH
      character (LEN=8) DAY
      character MISFIT_OPT
      logical ISLOG, ApplyConditions
      integer :: ierr, color, DA_Nprocs
      integer :: biol, bphy, nchl, sat, argo, SysErr, system
      real*8  :: chl_dep


      NAMELIST /biolst/ biol, bphy, nchl, chl_dep, sat, argo

      DAparttime  = MPI_WTIME()
      MONTH=datestr(5:6)
      DAY  =datestr(1:8)

      MISFIT_OPT ='2'
      ISLOG = .false.
      ApplyConditions = .false. ! .true.

      if(rank .eq. 0) then
        open(11,file='var_3d_nml',form='formatted')
        read(11,biolst)
        close(11)
      endif

      ! SATDIR li vuole pronti all'uso, già tagliati e interpolati


      SUFFIX   = '_d-OC_CNR-L4-CHL-MedOC4_SAM_7KM-MED-REP-v02.nc'
      SATFILE   = 'SATELLITE/' // DAY // trim(SUFFIX)
      VARFILE   = 'DA_static_data/MISFIT/VAR2D/var2D.' // MONTH // '.nc'
      EOF_FILE  = 'DA_static_data/3D_VAR/EOF/eof.'  // MONTH // '.nc'
      GRID_FILE = 'DA_static_data/3D_VAR/GRID/BFM_grid.nc'
      RCORR_FILE = 'DA_static_data/3D_VAR/gradsal.nc'

      MISFIT_FILE='DA__FREQ_1/'// DAY // '.chl_mis.nc'
      ARGO_FILE='DA__FREQ_1/'// DAY // '.P_l_arg_mis.dat'
      CORR_FILE = 'DA__FREQ_1/'// DAY // '_corr.nc'
      EIV_FILE  = 'DA__FREQ_1/'// DAY // '_eiv.nc'
      OBS_FILE = 'obs_1.dat' ! 'obs_'//fgrd//'.dat'


      CHLSUP_FOR_DA = 'DA__FREQ_1/chl.' // datestr // '.nc'

      CALL trcditDA(DATEstr, datefrom, DATEstr)! average of 12 h
      ! Sviluppo : se non serve avere il file ave delle DA, la snutell puo' fare il collecting
      CALL trcwriDA(DATEstr)  ! Dumps Before Assimilation real*4

      DA_Nprocs = 20
      if (rank .lt. DA_Nprocs ) then
          if(rank .eq. 0 .and. sat .eq. 1) then
            call CREATEMISFIT(SATFILE,VARFILE,MISFIT_OPT, ISLOG, MISFIT_FILE) ! produces MISFIT.nc
            write(*,*) 'eof = ',   trim(EOF_FILE)
            write(*,*) 'grid = ',  trim(GRID_FILE)
          endif

          ! if(rank .eq. 0) call system(ScriptName) !//" -t "//DAY)
          if(rank .eq. 0) then
            SysErr = system("../float_preproc/Float_misfit_gen.sh -d ../float_preproc -t "//DAY)
            if(SysErr /= 0) call MPI_Abort(MPI_COMM_WORLD, -1, SysErr)
          endif
  
          call MPI_Barrier(Var3DCommunicator, ierr)


          call OCEANVAR
          if(rank .eq. 0) CALL SNUTELL(datestr, ISLOG, ApplyConditions)
      endif

      call mppsync()
      CALL trcrstDA(datestr)

      DAparttime  = MPI_WTIME() - DAparttime
      DAtottime   = DAtottime + DAparttime

      ! Sviluppo : se non serve il RST After Assimilation, puo' fare uno scatter della nuova variabile di stato
      END SUBROUTINE mainASSIMILATION
