
!*********************************************************************
!                   Institut Pierre-Simon Laplace
!
!                             L O D Y C
!       Laboratoire d Oceanographie DYnamique et de CLImatologie
!                   Unite mixte de recherche 121
!                          CNRS/ORSTOM/UPCM
!
!
!                             O P A model
!
!                    Ocean General Circulation Model
!
!
! version OPA.8.1  (1997)
!
!*********************************************************************
!
!  opa                 main program
!  **first part : initialization
!          inipar              initialization parameters
!          parcst            constants
!           parlec            namelist read
!          parctl            control options
!                 ctlopn              open file and check
!     inimpp                          Initialisation for large domain
!     inimpp2                         Initialisation with land elimination
!     domrea                         read a meshmask file
!     domstp                         prepares diferent time step for
!                                     each level
!     initrc                            initialization passive tracers
!          ......            passive tracers model package
!     inihdf                          horizontal diffusion coefficient
!     dtadyn                          dynamique and physique
!    mlbxinit            task initialization
!     **second part : time step loop
!    step                 time step routine
!           dtadyn            dynamique and physique
!                 div                horizontal divergence
!                 wzv                vertical velocity
!                 eos            equation of state
!                 bn2            Brunt-Vaisala frequency
!                 hdfslp                isopycnal slope computation
!           stpctl            control step execution
!                 mpprisl             mpp real(8) no zero value
!           trcwri                     write restart FILE for tracer
!           trcdia            write standard output
!     **end of time loop
!
!*********************************************************************
!
      SUBROUTINE myopa()
!---------------------------------------------------------------------
!
!                        ******************
!                        *   ***********  *
!                        *    OPA MODEL   *
!                        *   ***********  *
!                        ******************
!
!  Purpose :
!  --------
!     opa solves the primitive equations with rigid lid assumption
!    on any orthogonal curvilinear mesh on the sphere which vertical
!    direction is the geographical one.
!
!   Method :
!   -------
!    model general initialization
!
!    execution of time loop by calling the routine step
!
!   Input :        : no
!   ------
!
!   External :
!   --------
!      inipar, inidta, inidtr,
!    step.

!
! parameters and commons
! ======================

       USE myalloc
       USE myalloc_mpp
       USE IO_mem
       USE FN_mem
       USE ADV_mem
       USE HDF_mem
       USE ZDF_mem
       USE OPT_mem
       USE BC_mem
       USE BIO_mem
       USE time_manager

!----------------------------------------------------------------------
! externals
! =========
!      EXTERNAL step
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER ME
      INTEGER mynode

! 0. Initializations
! ==================

      numout=2  ! 0.1 Initialization of unit print

      ! *********************************************

      OPEN(UNIT=numout,FILE='ocean.output',FORM='FORMATTED')
      ME = mynode() !  Nodes selection
      call parini(ME)

      narea = ME+1
      lwp = narea.EQ.1


      call parlec  ! read namelist
      call time_init

      write(*,*) 'time_init called'

      ! -------------------------
      call ALLOC_ALL ! Needs Time_Manager
      ! -------------------------

      numnam = 208
      numrst = 209
      numbat = 212
      nummsh = 214 ! serve alla restget
      nummpp =  8
       nimpp  = 1
       njmpp  = 1 

#ifdef key_mpp
      OPEN(UNIT=nummpp,FILE='mpp.output',FORM='FORMATTED')
#endif



! 0.3 First prints

      IF(lwp) THEN
          WRITE(numout,*) ' '
          WRITE(numout,*) '             Institut Pierre-Simon LAPLACE'
          WRITE(numout,*) '                      L O D Y C'
          WRITE(numout,*) ' Laboratoire d"Oceanographie DYnamique et de Climatologie'
          WRITE(numout,*) '             unite mixte de recherche 121'
          WRITE(numout,*) '                  CNRS/ORSTOM/UPCM'
          WRITE(numout,*) ' '
          WRITE(numout,*) ' '
          WRITE(numout,*) '                  Offline O P A model'
          WRITE(numout,*) ' '
          WRITE(numout,*) '            Ocean General Circulation Model'
          WRITE(numout,*) ' '
          WRITE(numout,*) ' '
          WRITE(numout,*) ' version OPA.8.1.offline'
          WRITE(numout,*) ' 03/04/1998'
          WRITE(numout,*) ' '
          WRITE(numout,*) ' '
      ENDIF




      CALL INITIALIZING

! -----------------------------------------------
      CALL step
! -----------------------------------------------


! 3. Job end
! ==========

#ifdef key_mpp
      write(*,*)'My_Rank=',Rank,'before mpi finalize'
      CALL mppstop
      write(*,*)'My_Rank=',Rank,'after mpi finalize'
#endif

      if(lwp) WRITE(numout,*) 'End of calculation. Good bye.'

      CLOSE( numout ) ! others units are closed in mppstop

      END SUBROUTINE myopa





! ***************************************************************
! ***************************************************************
      SUBROUTINE ALLOC_ALL

       USE myalloc
       USE myalloc_mpp
       USE IO_mem
       USE FN_mem
       USE ADV_mem
       USE HDF_mem
       USE ZDF_mem
       USE OPT_mem
       USE BC_mem
       USE BIO_mem
       USE SED_mem

       IMPLICIT NONE

       REAL(8) mem_all_tot
       INTEGER ierr

      write(*,*)'My_Rank=',Rank,':allocazione_memoria'
      mem_all_tot=0
      mem_all=0

      call   alloc_tot() ; write(*,*)'My_Rank:',Rank,'alloc_init:', mem_all ; mem_all_tot=mem_all_tot+mem_all
      call  myalloc_DR() ; write(*,*)'My_Rank:',Rank,'alloc_DR  :', mem_all ; mem_all_tot=mem_all_tot+mem_all
      call myalloc_OPT() ; write(*,*)'My_Rank:',Rank,'alloc_OPT :', mem_all ; mem_all_tot=mem_all_tot+mem_all
      call myalloc_ADV() ; write(*,*)'My_Rank:',Rank,'alloc_ADV :', mem_all ; mem_all_tot=mem_all_tot+mem_all
      call myalloc_HDF() ; write(*,*)'My_Rank:',Rank,'alloc_HDF :', mem_all ; mem_all_tot=mem_all_tot+mem_all
      call myalloc_ZDF() ; write(*,*)'My_Rank:',Rank,'alloc_ZDF :', mem_all ; mem_all_tot=mem_all_tot+mem_all
      call myalloc_DYN() ; write(*,*)'My_Rank:',Rank,'alloc_DYN :', mem_all ; mem_all_tot=mem_all_tot+mem_all

#ifdef key_trc_dmp
!     needs Time_Manager
      call alloc_DTATRC(); write(*,*)'My_Rank:',Rank,'alloc_TRC:', mem_all ; mem_all_tot=mem_all_tot+mem_all
#endif

      call myalloc_BIO() ; write(*,*)'My_Rank:',Rank,'alloc_BIO :', mem_all ; mem_all_tot=mem_all_tot+mem_all
      call myalloc_SED() ; write(*,*)'My_Rank:',Rank,'alloc_SED :', mem_all ; mem_all_tot=mem_all_tot+mem_all

      call myalloc_FN()  ; write(*,*)'My_Rank:',Rank,'alloc_FN  :', mem_all ; mem_all_tot=mem_all_tot+mem_all


      call MPI_ALLREDUCE(jpi, jpi_max, 1, MPI_INTEGER, MPI_MAX,MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(jpj, jpj_max, 1, MPI_INTEGER, MPI_MAX,MPI_COMM_WORLD, ierr)

      call myalloc_IO()  ; write(*,*)'My_Rank:',Rank,'alloc_IO  :', mem_all ; mem_all_tot=mem_all_tot+mem_all

      write(*,*)'My_Rank,',Rank,'memoria_totale_allocata',mem_all_tot
      END SUBROUTINE ALLOC_ALL



! *************************************************************
!      SUBROUTINE INITIALIZING
! *************************************************************
       SUBROUTINE INITIALIZING
       USE myalloc
       IMPLICIT NONE
! 1. Model general initialization
! ===============================

!    Run parameters
!    --------------
      CALL parcst;  CALL parlec ; CALL parctl ! CALL inipar
      CALL dialec


! ... Large domain
      IF ( jpni*jpnj.EQ.jpnij ) THEN
          CALL inimpp !   ... standard cutting out
        ELSE
!   ... eliminate land processors

!         CALL inimpp2 ! to be tested
         STOP
      ENDIF



      write(*,*) 'Startting domrea '
      CALL domrea   !   Domain
      CALL domstp   !   Time Step for each level for tracer
      CALL inihdf   !   diffusion horizontal coefficient

!    Passive tracers


      IF(lwp) THEN
          WRITE(numout,*) ' '
          WRITE(numout,*) ' *** number of passive tracer jptra = ',jptra
          WRITE(numout,*) ' '
      ENDIF

      CALL trclec ! READ passive tracers namelists
      CALL trcctl ! control consistency between parameters, cpp key and namelists
!      CALL time_init


      call trccof ! call trcdta ! initialisation of data fields
      CALL trcini ! computes some initialisation

      IF(lrsttr) THEN ! restart from a FILE (nutrst)
          CALL trcrst
      ELSE
          CALL trcdtr ! start from anything ELSE
      ENDIF



      call ave_init

      END SUBROUTINE INITIALIZING




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
       if (lwp) write(*,*) 'start End ok'
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


        call TimeInterpolation(sec,TC_FOR, TC_FOR%Before, TC_FOR%After, t_interp)
        call TimeInterpolation(sec,TC_TIN, TC_TIN%Before, TC_TIN%After, t_interp)
        call TimeInterpolation(sec,TC_ATM, TC_ATM%Before, TC_ATM%After, t_interp)
        call TimeInterpolation(sec,TC_GIB, TC_GIB%Before, TC_GIB%After, t_interp)
        call TimeInterpolation(sec,TC_LEX, TC_LEX%Before, TC_LEX%After, t_interp)

        if (lwp) then
            write(*,*) 'BeforeForcings', TC_FOR%Before, 'AfterForcing', TC_FOR%After
            write(*,*) 'BeforeRivers',   TC_TIN%Before, 'AfterRivers',  TC_TIN%After
            write(*,*) 'BeforeGib',      TC_GIB%Before, 'AfterGib',     TC_GIB%After
            write(*,*) 'BeforeAtm',      TC_ATM%Before, 'AfterAtm',     TC_GIB%After
            write(*,*) 'BeforeKex',      TC_LEX%Before, 'AfterKex',     TC_LEX%After
        endif
      END SUBROUTINE time_init


! *************************************************************
! ******** ave_init ******************************************
! *************************************************************

      SUBROUTINE ave_init

      USE myalloc
#ifdef key_mpp
      USE myalloc_mpp
#endif


      IMPLICIT NONE

          call mppsync()

          traIO     = trn ! initialization of Average  Matrix
          snIO      = sn
          tnIO      = tn
          vatmIO    = vatm
          empIO     = emp
          qsrIO     = qsr
          unIO      = un
          bblxIO    = bblx
          vnIO      = vn
          bblyIO    = bbly
          wnIO      = wn
          avtIO     = avt
          tra_ppIO = tra_pp


      end SUBROUTINE ave_init
