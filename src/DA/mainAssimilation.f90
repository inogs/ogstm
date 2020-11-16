      SUBROUTINE mainASSIMILATION(datestr,datefrom)

      use filenames
      use DA_mem
      use DA_Params
      use ogstm_mpi_module
      use mpi_str, only: Var3DCommunicator
      use calendar
      USE TREd_var_MP

      IMPLICIT NONE

      CHARACTER(LEN=17), INTENT(IN) :: datestr, dateFrom

      character(LEN=1024) SATFILE, VARFILE
      character(LEN=2) MONTH
      character (LEN=8) DAY
      character MISFIT_OPT
      logical ISLOG
      integer :: ierr, color, ii
      integer :: SysErr, system


      DAparttime  = MPI_WTIME()
      MONTH=datestr(5:6)
      DAY  =datestr(1:8)

      ! DA_Params init
      DA_Date = datestr
      ShortDate = datestr(1:11)//datestr(13:14)//datestr(16:17)
      jpk_200 = AssimilationLevels
      NBioVar = 17
      DA_JulianDate = datestring2sec(datestr)

      MISFIT_OPT ='2'
      ISLOG = .false.

      ! SATDIR li vuole pronti all'uso, già tagliati e interpolati
      SATFILE   = 'SATELLITE/' // DAY // trim(satfile_suffix)
      VARFILE   = 'DA_static_data/MISFIT/VAR2D/var2D.' // MONTH // '.nc'
      EOF_FILE  = 'DA_static_data/3D_VAR/EOF/eof.'  // MONTH // '.nc'
      GRID_FILE = 'DA_static_data/3D_VAR/GRID/BFM_grid.nc'
      ANIS_FILE = 'DA_static_data/3D_VAR/gradsal.nc'

      MISFIT_FILE='DA__FREQ_1/'// DAY // '.chl_mis.nc'
      ARGO_FILE='DA__FREQ_1/'// DAY // '.P_l_arg_mis.dat'
      CORR_FILE = 'DA__FREQ_1/'// DAY // '_corr.nc'
      EIV_FILE  = 'DA__FREQ_1/'// DAY // '_eiv.nc'
      OBS_FILE = 'obs_1.dat' ! 'obs_'//fgrd//'.dat'


      CHLSUP_FOR_DA = 'DA__FREQ_1/chl.' // datestr // '.nc'
      CALL trcwriDA(DATEstr)  ! Dumps Before Assimilation real*4

      if (V3D_VAR_PARALLEL) then

          allocate(DA_VarList(NBioVar))
          do ii=1,NBioVar
            DA_VarList(ii) = varlistDA(ii)
          enddo

          if(myrank .eq. 0) then ! .and. sat .eq. 1) then
            write(*,*) 'satfile=', trim(SATFILE)
            write(*,*) 'varfile=', trim(VARFILE)
            write(*,*) 'misfit=', trim(MISFIT_FILE)
            call CREATEMISFIT(SATFILE,VARFILE,MISFIT_OPT, ISLOG, MISFIT_FILE) ! produces MISFIT.nc
            write(*,*) 'eof = ',   trim(EOF_FILE)
            write(*,*) 'grid = ',  trim(GRID_FILE)
          endif

          ! if(myrank .eq. 0) then
          !   SysErr = system("../float_preproc/Float_misfit_gen.sh -d ../float_preproc -t "//DAY)
          !   if(SysErr /= 0) call MPI_Abort(MPI_COMM_WORLD, -1, SysErr)
          ! endif
  
          call MPI_Barrier(Var3DCommunicator, ierr)


          call BIOVAR

          ! deallocation of DA_VarList array
          call clean_da_params
      endif

      call mppsync()
      CALL trcrstDA(datestr)

      DAparttime  = MPI_WTIME() - DAparttime
      DAtottime   = DAtottime + DAparttime

      ! Sviluppo : se non serve il RST After Assimilation, puo' fare uno scatter della nuova variabile di stato
      END SUBROUTINE mainASSIMILATION
