! developed by Simone Spada (sspada@ogs.it) at OGS

module Ens_Mem
    USE mpi
    
    use myalloc, &
            only: lwp    
    use Time_Manager, &
        only: dump_container, &
            Load_Dump_container, Unload_Dump_container
    USE ogstm_mpi_module, &
            ONLY: glcomm
    
implicit none
    integer, parameter :: EnsRankZero=0, myrankZero=0, EnsIOUnit=86
    double precision, parameter :: Ens_Miss_val=1.0d20
    
    integer :: EnsDebug
    integer :: EnsComm, EnsRank, EnsSize
    logical :: UseParams
    integer :: EnsShareRestart, EnsShareParams
    logical :: EnsSaveEachRestart, EnsSaveMeanRestart, EnsSaveEachAve, EnsSaveMeanAve, EnsAveDouble
    character(Len=100) :: Ens_restart_prefix, Ens_restart_ens_prefix, &
        Ens_ave_freq_1_prefix, Ens_ave_freq_1_ens_prefix, Ens_ave_freq_2_prefix, Ens_ave_freq_2_ens_prefix, &
        Ens_flux_prefix, Ens_flux_ens_prefix

    type(dump_container) :: ForecastTimes, daTimes
    logical :: EnsSaveAfterForecast, EnsSaveAfterAnalysis
    integer :: LocalRange
    logical :: UseLocalObsDumping
    double precision :: ForgettingFactor
    character(len=100) :: Ens_forecast_prefix, Ens_forecast_ens_prefix, Ens_analysis_prefix, Ens_analysis_ens_prefix
    
    double precision, allocatable, dimension(:) :: EnsWeights
    
contains

subroutine Ens_Namelist(filename)
        
    character(len=*) :: filename
    
    NAMELIST/Ensemble_setup/ EnsDebug, EnsSize, UseParams, EnsShareRestart, EnsShareParams, &
        EnsSaveEachRestart, EnsSaveMeanRestart, EnsSaveEachAve, EnsSaveMeanAve, EnsAveDouble, &
        EnsSaveAfterForecast, EnsSaveAfterAnalysis, &
        Ens_restart_prefix, Ens_restart_ens_prefix, &
        Ens_ave_freq_1_prefix, Ens_ave_freq_1_ens_prefix, Ens_ave_freq_2_prefix, Ens_ave_freq_2_ens_prefix, &
        Ens_flux_prefix, Ens_flux_ens_prefix
        
    NAMELIST/Ensemble_DA_setup/ EnsSaveAfterForecast, EnsSaveAfterAnalysis, &
        LocalRange, UseLocalObsDumping, &
        ForgettingFactor, &
        Ens_forecast_prefix, Ens_forecast_ens_prefix, Ens_analysis_prefix, Ens_analysis_ens_prefix
    
    if (lwp) then
        
    end if
    
    OPEN(unit=EnsIOUnit, file=filename, status='OLD')
    
        EnsDebug=0
        EnsSize=1
        UseParams=.false.

        !*****************************************!
        ! This part is relevant only for EnsSize>1
        EnsShareRestart=0
        EnsShareParams=0
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

        REWIND(EnsIOUnit)
        READ(EnsIOUnit, Ensemble_setup)
        
        if (EnsSize<1) EnsSize=1

        IF(lwp) THEN
            write(*,*) 'Namelist Ensemble parameters:'
            WRITE(*,*) ' '
            WRITE(*,*) 'Ensemble_setup'
            WRITE(*,*) ' '
            WRITE(*,*) ' EnsDebug (verbosity): ', EnsDebug
            WRITE(*,*) ' EnsSize: ', EnsSize
            WRITE(*,*) ' UseParams: ', UseParams
            WRITE(*,*) ' EnsShareRestart: ', EnsShareRestart
            WRITE(*,*) ' EnsShareParams: ', EnsShareParams
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
        LocalRange=30
        UseLocalObsDumping=.true.
        ForgettingFactor=0.9d0
        Ens_forecast_prefix='ENS_FORECAST/RST'
        Ens_forecast_ens_prefix='ENS_FORECAST/ENSEMBLE/RST'
        Ens_analysis_prefix='ENS_ANALYSIS/RST'
        Ens_analysis_ens_prefix='ENS_ANALYSIS/ENSEMBLE/RST'

        REWIND(EnsIOUnit)
        READ(EnsIOUnit, Ensemble_DA_setup)

        IF(lwp) THEN
            WRITE(*,*) ' '
            WRITE(*,*) 'Ensemble_DA_setup'
            WRITE(*,*) ' '
            WRITE(*,*) ' EnsSaveAfterForecast: ', EnsSaveAfterForecast
            WRITE(*,*) ' EnsSaveAfterAnalysis: ', EnsSaveAfterAnalysis
            WRITE(*,*) ' LocalRange: ', LocalRange
            WRITE(*,*) ' UseLocalObsDumping: ', UseLocalObsDumping
            WRITE(*,*) ' ForgettingFactor: ', ForgettingFactor
            WRITE(*,*) ' Ens_forecast_prefix: ', TRIM(Ens_forecast_prefix)
            WRITE(*,*) ' Ens_forecast_ens_prefix: ', trim(Ens_forecast_ens_prefix)
            WRITE(*,*) ' Ens_analysis_prefix: ', trim(Ens_analysis_prefix)
            WRITE(*,*) ' Ens_analysis_ens_prefix: ', trim(Ens_analysis_ens_prefix)
        END IF

    CLOSE(EnsIOUnit)
    
end subroutine

subroutine Ens_allocate
    
    if (EnsSize>1) then
        
        allocate(EnsWeights(0:EnsSize-1))
        EnsWeights=1.0d0/EnsSize

#ifdef ExecEnsDA
        ForecastTimes%FileName = 'forecastTimes'
        ForecastTimes%Name='...'
        call Load_Dump_container(ForecastTimes)
        
        daTimes%FileName = 'daTimes'
        daTimes%Name='...'
        call Load_Dump_container(daTimes)
#endif
    end if

end subroutine

subroutine Ens_deallocate
    
    if (EnsSize>1) then
    
        deallocate(EnsWeights)
        
#ifdef ExecEnsDA
        call Unload_Dump_container(ForecastTimes)
        call Unload_Dump_container(daTimes)
#endif
    end if

end subroutine

Subroutine Ens_shared_alloc(member_size, member_pointer, global_pointer, window)        
    USE, INTRINSIC :: ISO_C_BINDING
    
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

    call MPI_Win_shared_query(window, MPI_PROC_NULL, window_buffer_size, double_size, C_pointer, ierror)
    
    CALL C_F_POINTER(C_pointer, global_pointer, SHAPE = (/member_size, EnsSize/))
end subroutine

Subroutine Ens_shared_alloc_base(member_size, member_pointer, global_pointer, window)
    USE, INTRINSIC :: ISO_C_BINDING
    
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
    
    if (EnsRank==EnsRankZero) then
        CALL MPI_Win_allocate_shared(0, double_size, MPI_INFO_NULL, EnsComm, C_pointer, window, ierror)
    else
        CALL MPI_Win_allocate_shared(window_buffer_size, double_size, MPI_INFO_NULL, EnsComm, C_pointer, window, ierror)
    end if
    if (ierror/=0) then
        WRITE(*, *) 'Error with MPI_Win_allocate_shared. Aborting.'
        call MPI_abort(glcomm, 1, ierror)
    end if

    if (EnsRank/=EnsRankZero) CALL C_F_POINTER(C_pointer, member_pointer, SHAPE = (/member_size/))

    call MPI_Win_shared_query(window, MPI_PROC_NULL, window_buffer_size, double_size, C_pointer, ierror)
    
    CALL C_F_POINTER(C_pointer, global_pointer, SHAPE = (/member_size, EnsSize-1/))
end subroutine

LOGICAL FUNCTION IsaForecast(datestring)
    use Time_Manager, &
        only: DATESTART
    
    CHARACTER(LEN=17), INTENT(IN) :: datestring
    ! LOCAL
    INTEGER I

    IsaForecast = .false.
    
    !if (datestring.eq.DATESTART) return
    
    DO I=1, ForecastTimes%N
        if (datestring.eq.ForecastTimes%TimeStrings(I)) then
            IsaForecast = .true.
            return
        endif
    ENDDO
    
end function       

LOGICAL FUNCTION IsaDataAssimilation(datestring)
    use Time_Manager, &
        only: DATESTART
    
    CHARACTER(LEN=17), INTENT(IN) :: datestring
    ! LOCAL
    INTEGER I

    IsaDataAssimilation = .false.
    
    !if (datestring.eq.DATESTART) return
    
    DO I=1, daTimes%N
        if (datestring.eq.daTimes%TimeStrings(I)) then
            IsaDataAssimilation = .true.
            return
        endif
    ENDDO
    
end function       
    
end module
