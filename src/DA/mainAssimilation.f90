      SUBROUTINE mainASSIMILATION(datestr,datefrom)

      use filenames
      use DA_mem
! ---- modules from 3Dvar
      use DA_Params
      use drv_str
! -----------------------
      use ogstm_mpi_module
      use mpi_str, only: Var3DCommunicator
      use calendar, only: datestring2sec
      USE TREd_var_MP


      IMPLICIT NONE

      CHARACTER(LEN=17), INTENT(IN) :: datestr, dateFrom

      character(LEN=1024) SATFILE, VARFILE
      character(LEN=1024) HOSTDIR
      character(LEN=2) MONTH,HOUR
      character (LEN=8) DAY
      character (LEN=5) OBStype
      character MISFIT_OPT
      logical ISLOG
      integer :: ierr, color, ii
      integer :: SysErr, system
      integer :: number_of_argo_obs
      
      logical existFile, condition_to_3dvar

      DAparttime  = MPI_WTIME()
      MONTH=datestr(5:6)
      DAY  =datestr(1:8)
      HOUR =datestr(10:11)

      ! DA_Params init
      DA_Date = datestr
      ShortDate = ConvertDate(datestr)
      if(HOUR.eq.'12') then
        jpk_200 = AssimilationLevels_sat
        OBStype = '__Sat'
      elseif(HOUR.eq.'13') then
        jpk_200 = AssimilationLevels_float
        OBStype = 'Float'
      endif
      !NBioVar = 17
      DA_JulianDate = datestring2sec(DA_Date)

      MISFIT_OPT ='2'
      ISLOG = .false.

      ! SATDIR li vuole pronti all'uso, già tagliati e interpolati
      SATFILE   = 'SATELLITE/' // DAY // trim(satfile_suffix)
      VARFILE   = 'DA_static_data/VAR_SAT/var2D.' // MONTH // '.nc'
      !EOF_FILE_CHL  = 'DA_static_data/3D_VAR/EOF/CHL/eof.'  // MONTH // '.nc'
      EOF_FILE_CHL  = 'DA_static_data/3D_VAR/EOF/CHL' // OBStype // '/eof.'  // MONTH // '.nc'
      EOF_FILE_N3N  = 'DA_static_data/3D_VAR/EOF/N3n/eof.'  // MONTH // '.nc'
      EOF_FILE_O2O  = 'DA_static_data/3D_VAR/EOF/O2o/eof.'  // MONTH // '.nc'
      EOF_FILE_MULTI= 'DA_static_data/3D_VAR/EOF/CHLN3n/eof.'  // MONTH // '.nc'
      STD_FILE_MULTI= 'DA_static_data/3D_VAR/STD/std.'  // MONTH // '.nc'
      NUTCOV_FILE   = 'DA_static_data/3D_VAR/CROSSCORRS/crosscorrs.' // MONTH // '.nc'
      NUTCHLCOV_FILE   = 'DA_static_data/3D_VAR/CROSSCORRS/crosscorrs.' // MONTH // '.nc'

      GRID_FILE = 'DA_static_data/3D_VAR/GRID/BFM_grid' //OBStype// '.nc'
      ANIS_FILE = 'DA_static_data/3D_VAR/gradsal.nc'
      RCORR_FILE = 'DA_static_data/3D_VAR/chl_rad_corr.nc'

      MISFIT_FILE='DA__FREQ_1/'// DAY // '.chl_mis.nc'
      ARGO_FILE='DA__FREQ_1/'// DAY // '.arg_mis.dat'
      CORR_FILE = 'DA__FREQ_1/'// DAY // '_corr.nc'
      EIV_FILE  = 'DA__FREQ_1/'// DAY // '_eiv.nc'
      OBS_FILE = 'obs_1.dat' ! 'obs_'//fgrd//'.dat'


      CHLSUP_FOR_DA = 'DA__FREQ_1/chl.' // datestr // '.nc'
      CALL trcwriDA(DATEstr)  ! Dumps Before Assimilation real*4

      if (V3D_VAR_PARALLEL) then

!          allocate(DA_VarList(NBioVar))
!          do ii=1,NBioVar
!            DA_VarList(ii) = varlistDA(ii)
!          enddo

          IF(myrank .eq. 0) then
            call def_nml
            call def_nml_multi
            if (drv%sat_obs.eq.1) then
               write(*,*) '--- Preparing satellite misfit ---'
               call CREATEMISFIT(SATFILE,VARFILE,MISFIT_OPT, ISLOG, MISFIT_FILE) ! produces MISFIT.nc
               write(*,*) 'eof = ',   trim(EOF_FILE_CHL)
               write(*,*) 'grid = ',  trim(GRID_FILE)
             endif


            if (drv%argo_obs.eq.1) then
               write(*,*) '--- Preparing float misfit ---'
               write(*,*)'deleting end.txt and writing start.txt'
               open(unit=456,file="end.txt")
               close(456,status='delete')
               open(unit=567,file='start.txt')
               write(567,'(g0)') DAY
               close(567)

               do while(.TRUE.)
                  INQUIRE(FILE='end.txt', EXIST=existFile)
                  if (existFile) then
                     exit
                  end if
                  call sleep(2)
               end do

              call countline("DA__FREQ_1/"//DAY//".arg_mis.dat",number_of_argo_obs)
              condition_to_3dvar = number_of_argo_obs.gt.1
              if (.not.condition_to_3dvar) then
                  write(*,*) 'No argo obs available: skipping OCEANVAR.'
                  write(*,*) 'You can remove RSTbefore files by doing:'
                  write(*,*) 'rm -f DA__FREQ_1/RSTbefore.'//datestr//'*'
                  write(*,*) 'rm DA__FREQ_1/RSTbefore.'//datestr(1:11)//datestr(13:14)//datestr(16:17)//'*'
              endif
            endif


          ENDIF

          call MPI_Barrier(Var3DCommunicator, ierr)
          call MPI_Bcast(condition_to_3dvar,1,MPI_LOGICAL,0, Var3DCommunicator, ierr)


          if (condition_to_3dvar) call OCEANVAR


          ! deallocation of DA_VarList array
          !call clean_da_params
      endif

      call mppsync()

      call MPI_Bcast(condition_to_3dvar,1,MPI_LOGICAL,0, MPI_COMM_WORLD, ierr)
      if (condition_to_3dvar)  CALL trcrstDA(datestr)


      DAparttime  = MPI_WTIME() - DAparttime
      DAtottime   = DAtottime + DAparttime

      ! Sviluppo : se non serve il RST After Assimilation, puo' fare uno scatter della nuova variabile di stato
      END SUBROUTINE mainASSIMILATION
