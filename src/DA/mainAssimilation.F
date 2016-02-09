
      SUBROUTINE mainASSIMILATION(datestr,datefrom)

      use filenames
      use DA_mem
      IMPLICIT NONE

      CHARACTER(LEN=17), INTENT(IN) :: datestr, dateFrom

      character(LEN=1024) SATFILE, VARFILE
      character(LEN=46) SUFFIX
      character(LEN=2) MONTH
      character (LEN=8) DAY
      character MISFIT_OPT
      logical ISLOG


      DAparttime  = MPI_WTIME()
      MONTH=datestr(5:6)
      DAY  =datestr(1:8)

      MISFIT_OPT ='2'
      ISLOG = .false.

      ! SATDIR li vuole pronti all'uso, già tagliati e interpolati


      SUFFIX   = '_d-OC_CNR-L4-CHL-MedOC4_SAM_7KM-MED-REP-v02.nc'
      SATFILE   = 'SATELLITE/' // DAY // trim(SUFFIX)
      VARFILE   = 'DA_static_data/MISFIT/VAR2D/var2D.' // MONTH // '.nc'
      EOF_FILE  = 'DA_static_data/3D_VAR/EOF/eof.'  // MONTH // '.nc'
      GRID_FILE = 'DA_static_data/3D_VAR/GRID/BFM_grid.nc'

      MISFIT_FILE='DA__FREQ_1/'// DAY // '.chl_mis.nc'
      CORR_FILE = 'DA__FREQ_1/'// DAY // '_corr.nc'
      EIV_FILE  = 'DA__FREQ_1/'// DAY // '_eiv.nc'
      OBS_FILE = 'obs_1.dat' ! 'obs_'//fgrd//'.dat'


      CHLSUP_FOR_DA = 'DA__FREQ_1/chl.' // datestr // '.nc'

      CALL trcditDA(DATEstr, datefrom, DATEstr)! average of 12 h
      ! Sviluppo : se non serve avere il file ave delle DA, la snutell puo' fare il collecting
      CALL trcwriDA(DATEstr)  ! Dumps Before Assimilation real*4

      if (rank == 0 ) then
          call CREATEMISFIT(SATFILE,VARFILE,MISFIT_OPT, ISLOG, MISFIT_FILE) ! produces MISFIT.nc

          write(*,*) 'eof = ',   trim(EOF_FILE)
          write(*,*) 'grid = ',  trim(GRID_FILE)

          call OCEANVAR
          CALL SNUTELL(datestr, ISLOG)
      endif

      call mppsync()
      CALL trcrstDA(datestr)

      DAparttime  = MPI_WTIME() - DAparttime
      DAtottime   = DAtottime + DAparttime

      ! Sviluppo : se non serve il RST After Assimilation, puo' fare uno scatter della nuova variabile di stato
      END SUBROUTINE mainASSIMILATION
