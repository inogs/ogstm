
!*********************************************************************
!               Istituto Nazionale di Oceanografia e di
!                       Geofisica Sperimentale
!
!                                OGS
!
!
!
!                             OGSTM model
!
!
!     Giorgio Bolzon (gbolzon@ogs.trieste.it)
!     Paolo Lazzari  (plazzari@ogs.trieste.it)
!
!
!*********************************************************************
!

MODULE OGSTM

      USE myalloc
      ! epascolo USE myalloc_mpp
      USE IO_mem
      USE FN_mem
      USE ADV_mem
      !USE HDF_mem
      USE ZDF_mem
      USE OPT_mem
      USE BC_mem
      USE BIO_mem
      USE SED_mem
      USE DIA_mem
      USE CALENDAR
      USE time_manager
      USE mpi
      USE NODE_NAME
      USE NODES_MODULE
      USE MATRIX_VARS
      USE MPI_GATHER_INFO
      USE dtype_procs_string_module
      use module_step
      use api_bfm 
      USE TREd_var_MP


#ifdef Mem_Monitor
      USE check_mem
      USE iso_c_binding
#endif
#ifdef ExecDA
      USE DA_MEM
      USE DA_VARS_module
#endif

! ----------------------------------------------------------------------
!  BEGIN BC_REFACTORING SECTION
!  ---------------------------------------------------------------------

      use bc_set_mod

! ----------------------------------------------------------------------
!  END BC_REFACTORING SECTION
!  ---------------------------------------------------------------------

IMPLICIT NONE

CONTAINS

SUBROUTINE ogstm_launcher()

#ifdef key_trc_bfm
#include "BFM_module_list.h"
#endif

     
      double precision :: timetosolution

      
      timetosolution = MPI_Wtime()

      call ogstm_initialize()
      

      call step


      call ogstm_finalize()

      timetosolution = MPI_Wtime() - timetosolution
      print*,"TIME TO SOLUTION =",timetosolution

     
      END SUBROUTINE ogstm_launcher


! *************************************************************
!      SUBROUTINE ogstm_initialize
! *************************************************************
SUBROUTINE ogstm_initialize()

! local declarations
! ==================
     
      ! *********************************************

      CALL mynode() !  Nodes selection

      narea = myrank+1
      lwp = narea.EQ.1

      IF(lwp) THEN
          OPEN(UNIT=numout,FILE='ogstm.out',FORM='FORMATTED')
          WRITE(numout,*) ' '
          WRITE(numout,*) '          Istituto Nazionale di Oceanografia e di '
          WRITE(numout,*) '                  Geofisica Sperimentale'
          WRITE(numout,*) ' '
          WRITE(numout,*) '                           OGS'
          WRITE(numout,*) ' '
          WRITE(numout,*) ' '
          WRITE(numout,*) ' '
          WRITE(numout,*) '                        OGSTM model'
          WRITE(numout,*) ' '
          WRITE(numout,*) ' '
          WRITE(numout,*) ' '
          WRITE(numout,*) ' '
          WRITE(numout,*) '  Giorgio Bolzon (gbolzon@ogs.trieste.it)'
          WRITE(numout,*) '  Paolo Lazzari  (plazzari@ogs.trieste.it)'
          WRITE(numout,*) ' '
          WRITE(numout,*) ' '
      ENDIF

      call parini

      call parlec      ! read namelist.init
      call time_init
      call trclec

      ! -------------------------
      call ALLOC_ALL ! Needs Time_Manager
      ! -------------------------

!    Run parameters
!    --------------
      call parcst
      call parctl ! controls consistency between parameters, cpp key and namelists

! 1. Model general initialization
! ===============================

      call node_name_fill    ! give the name of each node

      call nodes_module_find    !calculate the number of nodes used

      call populate_matrix_vars   !define matrix of variables to dump

      call init_mpi_gather_info   ! initialise all processors inidices

      call DEFINE_3D_PARALLEL_PROCS ! initialise      

      call domrea        !   Domain

      call inihdf        !   diffusion horizontal coefficient

      call trccof        ! initialisation of data fields

      call set_to_zero() ! set to zero some data arrays

      call trcrst        ! read restarts

      call photo_init

#ifdef ExecDA
      call DA_Init

      call DA_VARS
#endif

      call init_phys

! Initialization of Biogeochemical reactor with 1D approach
      call BFM0D_NO_BOXES(jpk*jpj*jpi,jpi,jpj,jpk,jpi*jpj)
      parallel_rank=myrank
      call Init_bfm()
      call BFM0D_INIT_IO_CHANNELS()

      call Initialize()

END SUBROUTINE ogstm_initialize

! ***************************************************************
! ***************************************************************
SUBROUTINE ALLOC_ALL

       double precision  mem_all_tot
       INTEGER err, ierr

      mem_all_tot=0.
      mem_all=0

#ifdef Mem_Monitor

       mem_all = get_mem(err)
       print *,"mem_all",mem_all
#endif

       !write(*,*)'My_Rank=',myrank,': Memory Allocation - Basal - (MB):',  mem_all
       mem_all_tot=mem_all_tot+mem_all

      call   alloc_tot() 
       !write(*,*)'My_Rank:',myrank,'alloc_init (MB):', mem_all
       mem_all_tot=mem_all_tot+mem_all
      call myalloc_OPT() 
       !write(*,*)'My_Rank:',myrank,'alloc_OPT  (MB):', mem_all
       mem_all_tot=mem_all_tot+mem_all
      call myalloc_ADV() 
       !write(*,*)'My_Rank:',myrank,'alloc_ADV  (MB):', mem_all
       mem_all_tot=mem_all_tot+mem_all
      !call myalloc_HDF() 
      ! write(*,*)'My_Rank:',myrank,'alloc_HDF  (MB):', mem_all 
       mem_all_tot=mem_all_tot+mem_all
      call myalloc_ZDF() 
       !write(*,*)'My_Rank:',myrank,'alloc_ZDF  (MB):', mem_all
       mem_all_tot=mem_all_tot+mem_all


#ifdef key_trc_dmp
!     needs Time_Manager
      call alloc_DTATRC()
       !write(*,*)'My_Rank:',myrank,'alloc_TRC  (MB):', mem_all
       mem_all_tot=mem_all_tot+mem_all
#endif
      call alloc_DIA()   
       !write(*,*)'My_Rank:',myrank,'alloc_DIA  (MB):', mem_all
       mem_all_tot=mem_all_tot+mem_all

      call myalloc_BIO() 
       !write(*,*)'My_Rank:',myrank,'alloc_BIO  (MB):', mem_all
       mem_all_tot=mem_all_tot+mem_all
      call myalloc_SED() 
       !write(*,*)'My_Rank:',myrank,'alloc_SED  (MB):', mem_all
       mem_all_tot=mem_all_tot+mem_all

      call myalloc_FN()  
       !write(*,*)'My_Rank:',myrank,'alloc_FN   (MB):', mem_all
       mem_all_tot=mem_all_tot+mem_all
      
      
      call MPI_ALLREDUCE(jpi, jpi_max, 1, MPI_INTEGER, MPI_MAX,MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(jpj, jpj_max, 1, MPI_INTEGER, MPI_MAX,MPI_COMM_WORLD, ierr)

      call myalloc_IO()  
      ! LEVEL1 write(*,*)'My_Rank:',myrank,'alloc_IO   (MB):', mem_all 
      mem_all_tot=mem_all_tot+mem_all

      ! LEVEL1 write(*,*)'My_Rank,',myrank,'Total Allocated Memory (MB):',mem_all_tot

END SUBROUTINE ALLOC_ALL


! *************************************************************
! ******** time_init ******************************************
! *************************************************************

SUBROUTINE time_init

      double precision sec, t_interp

      DELTAT = rdt ! importing namelist value

!      call ioconf_calendar('gregorian')
      

      SELECT CASE (calendarType)
        CASE ( 1) 
       CALL ioconf_calendar('gregorian')
        CASE ( 0) 
       CALL ioconf_calendar('noleap')
        CASE (30) 
       CALL ioconf_calendar('360d')
      END SELECT




! -----------------------------------------
      call Load_Timestrings
      if (CheckStartEnd()) then
       if (lwp) write(*,*) 'start End ok ', DATESTART, ' ', DATE__END
        else
           if (lwp) write(*,*) 'Problems with start End. Program will stop. '
           STOP
      endif
! -----------------------------------------


        call getTimesteps(TimeStepStart, TimeStep__End)

        if (lwp) then
            write(*,*) 'Time Step Start : ', TimeStepStart
            write(*,*) 'Time Step End   : ', TimeStep__End
        endif


        sec        = TimeStepStart*deltaT + TIME_0


        call TimeExtension(DATESTART,TC_FOR)
        ! call TimeExtension(DATESTART,TC_TIN)
        call TimeExtension(DATESTART,TC_ATM)
        ! call TimeExtension(DATESTART,TC_GIB)
        call TimeExtension(DATESTART,TC_LEX)
        call TimeExtension(DATESTART,TC_CO2)


        call TimeInterpolation(sec,TC_FOR, TC_FOR%Before, TC_FOR%After, t_interp)
        ! call TimeInterpolation(sec,TC_TIN, TC_TIN%Before, TC_TIN%After, t_interp)
        call TimeInterpolation(sec,TC_ATM, TC_ATM%Before, TC_ATM%After, t_interp)
        ! call TimeInterpolation(sec,TC_GIB, TC_GIB%Before, TC_GIB%After, t_interp)
        call TimeInterpolation(sec,TC_LEX, TC_LEX%Before, TC_LEX%After, t_interp)
        call TimeInterpolation(sec,TC_CO2, TC_CO2%Before, TC_CO2%After, t_interp)

        if (lwp) then
            write(*,*) 'BeforeForcings', TC_FOR%Before, 'AfterForcing', TC_FOR%After
            ! write(*,*) 'BeforeRivers',   TC_TIN%Before, 'AfterRivers',  TC_TIN%After
            ! write(*,*) 'BeforeGib',      TC_GIB%Before, 'AfterGib',     TC_GIB%After
            write(*,*) 'BeforeAtm',      TC_ATM%Before, 'AfterAtm',     TC_ATM%After
            write(*,*) 'BeforeCo2',      TC_CO2%Before, 'AfterCo2',     TC_CO2%After
            write(*,*) 'BeforeKex',      TC_LEX%Before, 'AfterKex',     TC_LEX%After

        endif
END SUBROUTINE time_init

! *************************************************************

SUBROUTINE photo_init
       
       INTEGER :: ji,jj, julianday


      call tau2julianday(TimeStepStart, deltaT, julianday)
         do ji=1, jpi
      do jj =1, jpj
            DAY_LENGTH(jj,ji) = photoperiod(julianday, gphit(jj,ji))
         enddo
      enddo


END SUBROUTINE photo_init

SUBROUTINE set_to_zero()

! Physical arrays set to zero

      un        = 0.0
      vn        = 0.0
      wn        = 0.0
      avt       = 0.0
      tn        = 0.0
      sn        = 0.0
! Passive tracers arrays set to zero

      xpar      = 0.0
      trn       = 0.0
      tra       = 0.0

END SUBROUTINE set_to_zero

! *************************************************************
!      SUBROUTINE ogstm_finalize
! *************************************************************

SUBROUTINE ogstm_finalize()

      CALL mppstop

      if(lwp) WRITE(numout,*) 'End of calculation. Good bye.'

! ----------------------------------------------------------------------
!  BEGIN BC_REFACTORING SECTION
!  ---------------------------------------------------------------------

      call boundaries%bc_set_destructor()
      deallocate(boundaries)

! ----------------------------------------------------------------------
!  END BC_REFACTORING SECTION
!  ---------------------------------------------------------------------

      if (lwp) CLOSE( numout ) ! others units are closed in mppstop
      CLOSE( numnam )

      ! clean memory
      call clean_memory()
      call clean_memory_bio()
      call clean_memory_fn()
      call clean_memory_opt()
      call clean_memory_sed()
      call unload_timestrings()
      call clean_memory_bc()
      call clean_memory_dia()
      call clean_memory_io()
      call clean_memory_adv()
      call clean_memory_zdf()

      END SUBROUTINE ogstm_finalize



END MODULE OGSTM

