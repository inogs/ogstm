! developed by Simone Spada (sspada@ogs.it) at OGS

module Ens_Mem
    use Time_Manager, &
        only: dump_container, &
            Load_Dump_container, Unload_Dump_container
    
implicit none
    integer, parameter :: EnsRankZero=0, EnsIOUnit=86
    double precision, parameter :: Ens_Miss_val=1.0d20
    
    integer :: EnsDebug, EnsShareRestart
    integer :: EnsComm, EnsRank, EnsSize
    logical :: EnsSaveEachRestart, EnsSaveMeanRestart, EnsSaveEachAve, EnsSaveMeanAve, EnsAveDouble
    character(Len=100) :: Ens_restart_prefix, Ens_restart_ens_prefix, &
        Ens_ave_freq_1_prefix, Ens_ave_freq_1_ens_prefix, Ens_ave_freq_2_prefix, Ens_ave_freq_2_ens_prefix, &
        Ens_flux_prefix, Ens_flux_ens_prefix
       
    logical :: EnsSaveAfterForecast, EnsSaveAfterAnalysis
    type(dump_container) :: ForecastTimes
    
    
contains

    subroutine Ens_Namelist(filename)
        use myalloc, &
            only: lwp
            
        character(len=*) :: filename
        
        NAMELIST/Ensemble_setup/ EnsDebug, EnsSize, EnsShareRestart, &
            EnsSaveEachRestart, EnsSaveMeanRestart, EnsSaveEachAve, EnsSaveMeanAve, EnsAveDouble, &
            EnsSaveAfterForecast, EnsSaveAfterAnalysis, &
            Ens_restart_prefix, Ens_restart_ens_prefix, &
            Ens_ave_freq_1_prefix, Ens_ave_freq_1_ens_prefix, Ens_ave_freq_2_prefix, Ens_ave_freq_2_ens_prefix, &
            Ens_flux_prefix, Ens_flux_ens_prefix
            
        NAMELIST/Ensemble_DA_setup/ EnsSaveAfterForecast, EnsSaveAfterAnalysis
        
        if (lwp) then
            
        end if
        
        OPEN(unit=EnsIOUnit, file=filename, status='OLD')
        
            EnsDebug=0
            EnsSize=1

            !*****************************************!
            ! This part is relevant only for EnsSize>1
            EnsShareRestart=1
            EnsSaveEachRestart=.true.
            EnsSaveMeanRestart=.true.
            EnsSaveEachAve=.true.
            EnsSaveMeanAve=.true.
            !*****************************************!

            EnsAveDouble=.false. ! saves ave faster and in double precision. Double space needed.
            Ens_restart_prefix='RESTART/RST'
            Ens_restart_ens_prefix='RESTART/ENSEMBLE/RST'
            Ens_ave_freq_1_prefix='AVE_FREQ_1/ave'
            Ens_ave_freq_1_ens_prefix='AVE_FREQ_1/ENSEMBLE/ave'
            Ens_ave_freq_2_prefix='AVE_FREQ_2/ave'
            Ens_ave_freq_2_ens_prefix='AVE_FREQ_2/ENSEMBLE/ave'
            Ens_flux_prefix='FLUXES/flux'
            Ens_flux_ens_prefix='FLUXES/ENSEMBLE/flux'

            !REWIND(EnsIOUnit)
            READ(EnsIOUnit, Ensemble_setup)
            
            if (EnsSize<1) EnsSize=1

            IF(lwp) THEN
                write(*,*) 'Namelist Ensemble parameters:'
                WRITE(*,*) ' '
                WRITE(*,*) 'Ensemble_setup'
                WRITE(*,*) ' '
                WRITE(*,*) ' EnsDebug (verbosity): ', EnsDebug
                WRITE(*,*) ' EnsSize: ', EnsSize
                WRITE(*,*) ' EnsShareRestart: ', EnsShareRestart
                WRITE(*,*) ' EnsSaveEachRestart: ', EnsSaveEachRestart
                WRITE(*,*) ' EnsSaveMeanRestart: ', EnsSaveMeanRestart
                WRITE(*,*) ' EnsSaveEachAve: ', EnsSaveEachAve
                WRITE(*,*) ' EnsSaveMeanAve: ', EnsSaveMeanAve
                WRITE(*,*) ' EnsAveDouble: ', EnsAveDouble
                WRITE(*,*) ' Ens_restart_prefix: ', trim(Ens_restart_prefix)
                WRITE(*,*) ' Ens_restart_ens_prefix: ', trim(Ens_restart_ens_prefix)
                WRITE(*,*) ' Ens_ave_freq_1_prefix: ', trim(Ens_ave_freq_1_prefix)
                WRITE(*,*) ' Ens_ave_freq_1_ens_prefix: ', trim(Ens_ave_freq_1_ens_prefix)
                WRITE(*,*) ' Ens_ave_freq_2_prefix: ', trim(Ens_ave_freq_2_prefix)
                WRITE(*,*) ' Ens_ave_freq_2_ens_prefix: ', trim(Ens_ave_freq_2_ens_prefix)
                WRITE(*,*) ' Ens_flux_prefix: ', trim(Ens_flux_prefix)
                WRITE(*,*) ' Ens_flux_ens_prefix: ', trim(Ens_flux_ens_prefix)
            END IF

            EnsSaveAfterForecast=.false.
            EnsSaveAfterAnalysis=.true.

            REWIND(EnsIOUnit)
            READ(EnsIOUnit, Ensemble_DA_setup)

            IF(lwp) THEN
                WRITE(*,*) ' '
                WRITE(*,*) 'Ensemble_DA_setup'
                WRITE(*,*) ' '
                WRITE(*,*) ' EnsSaveAfterForecast: ', EnsSaveAfterForecast
                WRITE(*,*) ' EnsSaveAfterAnalysis: ', EnsSaveAfterAnalysis
            END IF

        CLOSE(EnsIOUnit)
        
    end subroutine
    
    subroutine Ens_allocate
    
#ifdef ExecEnsDA 
        ForecastTimes%FileName = 'forecastTimes'
        ForecastTimes%Name='...'
        call Load_Dump_container(ForecastTimes)

#endif
    end subroutine
    
    subroutine Ens_deallocate
    
#ifdef ExecEnsDA 
        call Unload_Dump_container(ForecastTimes)

#endif
    end subroutine
    
    Subroutine Ens_shared_alloc(member_size, member_pointer, global_pointer, window)
        USE mpi
        USE, INTRINSIC :: ISO_C_BINDING
        
        USE ogstm_mpi_module, &
            ONLY: glcomm
        
        INTEGER, intent(in) :: member_size
        double precision, dimension(:), POINTER, contiguous, intent(out) :: member_pointer
        double precision, dimension(:,:), POINTER, contiguous, intent(out) :: global_pointer
        INTEGER, intent(out) :: window
        
        INTEGER :: ierror
        INTEGER :: double_size
        INTEGER(KIND=MPI_ADDRESS_KIND) :: window_buffer_size
        TYPE(C_PTR) :: C_pointer
        
        CALL MPI_Type_size(MPI_real8, double_size, ierror)
        window_buffer_size = double_size * member_size
    
        CALL MPI_Win_allocate_shared(window_buffer_size, double_size, MPI_INFO_NULL, EnsComm, C_pointer, window, ierror)
        if (ierror/=0) then
            WRITE(*, *) 'Error with MPI_Win_allocate_shared. Aborting.'
            call MPI_abort(glcomm, 1, ierror)
        end if
 
        CALL C_F_POINTER(C_pointer, member_pointer, SHAPE = (/member_size/))
    
        call MPI_Win_shared_query(window, 0, window_buffer_size, double_size, C_pointer, ierror)
        
        CALL C_F_POINTER(C_pointer, global_pointer, SHAPE = (/member_size, EnsSize/))
    end subroutine

    
    LOGICAL FUNCTION IsaForecast(datestring)
        use Time_Manager, &
            only: DATESTART
        
        CHARACTER(LEN=17), INTENT(IN) :: datestring
        ! LOCAL
        INTEGER I

        IsaForecast = .false.
        
        if (datestring.eq.DATESTART) return
        
        DO I=1, ForecastTimes%N
            if (datestring.eq.ForecastTimes%TimeStrings(I)) then
                IsaForecast = .true.
                return
            endif
        ENDDO
        
    end function
        
        
end module
