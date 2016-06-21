
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
      SUBROUTINE ogstm()
!



#include "OGSTM_module_list.h"
#ifdef key_trc_bfm
#include "BFM_module_list.h"
#endif



      call ogstm_initialize()


      call step


      call ogstm_finalize()



      END SUBROUTINE ogstm


! *************************************************************
!      SUBROUTINE ogstm_initialize
! *************************************************************
       SUBROUTINE ogstm_initialize()

#include "OGSTM_module_list.h"

       IMPLICIT NONE
! local declarations
! ==================
      INTEGER ME
      INTEGER mynode
      ! *********************************************

      OPEN(UNIT=numout,FILE='ocean.output',FORM='FORMATTED')

      ME = mynode() !  Nodes selection

      narea = ME+1
      lwp = narea.EQ.1

      IF(lwp) THEN
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

      call parini(ME)

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

      call inimpp

      call domrea        !   Domain

      call inihdf        !   diffusion horizontal coefficient

      call trccof        ! initialisation of data fields

      call set_to_zero() ! set to zero some data arrays

      call trcrst        ! read restarts

      call photo_init

! Initialization of Biogeochemical recator with 0D approach

      call BFM0D_NO_BOXES(1,1,1,1,1)
!     call BFM0D_INIT_IO_CHANNELS()
      call Initialize()

      END SUBROUTINE ogstm_initialize

! ***************************************************************
! ***************************************************************
      SUBROUTINE ALLOC_ALL

#include "OGSTM_module_list.h"

       IMPLICIT NONE

       REAL(8)  mem_all_tot
       INTEGER err, ierr

      mem_all_tot=0.
      mem_all=0

#ifdef Mem_Monitor
       mem_all = get_mem(err)
#endif

      write(*,*)'My_Rank=',Rank,': Memory Allocation - Basal - (MB):',  mem_all ; mem_all_tot=mem_all_tot+mem_all

      call   alloc_tot() ; write(*,*)'My_Rank:',Rank,'alloc_init (MB):', mem_all ; mem_all_tot=mem_all_tot+mem_all
      call myalloc_OPT() ; write(*,*)'My_Rank:',Rank,'alloc_OPT  (MB):', mem_all ; mem_all_tot=mem_all_tot+mem_all
      call myalloc_ADV() ; write(*,*)'My_Rank:',Rank,'alloc_ADV  (MB):', mem_all ; mem_all_tot=mem_all_tot+mem_all
      call myalloc_HDF() ; write(*,*)'My_Rank:',Rank,'alloc_HDF  (MB):', mem_all ; mem_all_tot=mem_all_tot+mem_all
      call myalloc_ZDF() ; write(*,*)'My_Rank:',Rank,'alloc_ZDF  (MB):', mem_all ; mem_all_tot=mem_all_tot+mem_all


#ifdef key_trc_dmp
!     needs Time_Manager
      call alloc_DTATRC(); write(*,*)'My_Rank:',Rank,'alloc_TRC  (MB):', mem_all ; mem_all_tot=mem_all_tot+mem_all
#endif
      call alloc_DIA()   ; write(*,*)'My_Rank:',Rank,'alloc_DIA  (MB):', mem_all ; mem_all_tot=mem_all_tot+mem_all

      call myalloc_BIO() ; write(*,*)'My_Rank:',Rank,'alloc_BIO  (MB):', mem_all ; mem_all_tot=mem_all_tot+mem_all
      call myalloc_SED() ; write(*,*)'My_Rank:',Rank,'alloc_SED  (MB):', mem_all ; mem_all_tot=mem_all_tot+mem_all

      call myalloc_FN()  ; write(*,*)'My_Rank:',Rank,'alloc_FN   (MB):', mem_all ; mem_all_tot=mem_all_tot+mem_all


      call MPI_ALLREDUCE(jpi, jpi_max, 1, MPI_INTEGER, MPI_MAX,MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(jpj, jpj_max, 1, MPI_INTEGER, MPI_MAX,MPI_COMM_WORLD, ierr)

      call myalloc_IO()  ; write(*,*)'My_Rank:',Rank,'alloc_IO   (MB):', mem_all ; mem_all_tot=mem_all_tot+mem_all

      write(*,*)'My_Rank,',Rank,'Total Allocated Memory (MB):',mem_all_tot
      END SUBROUTINE ALLOC_ALL


! *************************************************************
! ******** time_init ******************************************
! *************************************************************

      SUBROUTINE time_init

      USE calendar
      USE time_manager
      USE myalloc, ONLY: lwp, calendarType, rdt

      IMPLICIT NONE
      real(8) sec, t_interp

      DELTAT = rdt ! importing namelist value

!      call ioconf_calendar('gregorian');

      SELECT CASE (calendarType)
        CASE ( 1) ; CALL ioconf_calendar('gregorian')
        CASE ( 0) ; CALL ioconf_calendar('noleap')
        CASE (30) ; CALL ioconf_calendar('360d')
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
        call TimeExtension(DATESTART,TC_TIN)
        call TimeExtension(DATESTART,TC_ATM)
        call TimeExtension(DATESTART,TC_GIB)
        call TimeExtension(DATESTART,TC_LEX)
        call TimeExtension(DATESTART,TC_CO2)


        call TimeInterpolation(sec,TC_FOR, TC_FOR%Before, TC_FOR%After, t_interp)
        call TimeInterpolation(sec,TC_TIN, TC_TIN%Before, TC_TIN%After, t_interp)
        call TimeInterpolation(sec,TC_ATM, TC_ATM%Before, TC_ATM%After, t_interp)
        call TimeInterpolation(sec,TC_GIB, TC_GIB%Before, TC_GIB%After, t_interp)
        call TimeInterpolation(sec,TC_LEX, TC_LEX%Before, TC_LEX%After, t_interp)
        call TimeInterpolation(sec,TC_CO2, TC_CO2%Before, TC_CO2%After, t_interp)

        if (lwp) then
            write(*,*) 'BeforeForcings', TC_FOR%Before, 'AfterForcing', TC_FOR%After
            write(*,*) 'BeforeRivers',   TC_TIN%Before, 'AfterRivers',  TC_TIN%After
            write(*,*) 'BeforeGib',      TC_GIB%Before, 'AfterGib',     TC_GIB%After
            write(*,*) 'BeforeAtm',      TC_ATM%Before, 'AfterAtm',     TC_GIB%After
            write(*,*) 'BeforeCo2',      TC_CO2%Before, 'AfterCo2',     TC_CO2%After
            write(*,*) 'BeforeKex',      TC_LEX%Before, 'AfterKex',     TC_LEX%After

        endif
      END SUBROUTINE time_init



      SUBROUTINE photo_init
       USE myalloc
       USE calendar
       USE TIME_Manager
       IMPLICIT NONE
      ! LOCAL
       INTEGER ji,jj, julianday


      call tau2julianday(TimeStepStart, deltaT, julianday)
      do jj =1, jpj
         do ji=1, jpi
            DAY_LENGTH(ji,jj) = photoperiod(julianday, gphit(ji,jj))
         enddo
      enddo


      END SUBROUTINE photo_init

      SUBROUTINE set_to_zero()
#include "OGSTM_module_list.h"

       IMPLICIT NONE
! Physical arrays set to zero

      un        = 0.0
      vn        = 0.0
      wn        = 0.0
      avt       = 0.0
      tn        = 0.0
      sn        = 0.0
      hdivn     = 0.0
! Passive tracers arrays set to zero

      xpar      = 0.0
      trn       = 0.0
      tra       = 0.0

      END SUBROUTINE set_to_zero

! *************************************************************
!      SUBROUTINE ogstm_finalize
! *************************************************************
       SUBROUTINE ogstm_finalize()

#include "OGSTM_module_list.h"

       IMPLICIT NONE

      CALL mppstop

      if(lwp) WRITE(numout,*) 'End of calculation. Good bye.'

      CLOSE( numout ) ! others units are closed in mppstop
      CLOSE( numnam )

      END SUBROUTINE ogstm_finalize


