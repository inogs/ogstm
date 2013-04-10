
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
       USE DIA_mem
       USE CALENDAR
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

      narea = ME+1
      lwp = narea.EQ.1
      call parini(ME)

      call parlec  ! read namelist
      call time_init

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






      CALL INITIALIZING
      if(lwp) WRITE(numout,*) 'End of initialization. Good bye.'

      CLOSE( numout ) ! others units are closed in mppstop

! -----------------------------------------------
      CALL step
! -----------------------------------------------


! 3. Job end
! ==========

#ifdef key_mpp
      CALL mppstop
      if(lwp) write(*,*)'My_Rank=',Rank,'after mpi finalize'
#endif



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
       USE DIA_mem

       IMPLICIT NONE

       REAL(8) mem_all_tot
       INTEGER ierr

      mem_all_tot=0
      mem_all=0

      call   alloc_tot() ; mem_all_tot=mem_all_tot+mem_all
      call  myalloc_DR() ; mem_all_tot=mem_all_tot+mem_all
      call myalloc_OPT() ; mem_all_tot=mem_all_tot+mem_all
      call myalloc_ADV() ; mem_all_tot=mem_all_tot+mem_all
      call myalloc_HDF() ; mem_all_tot=mem_all_tot+mem_all
      call myalloc_ZDF() ; mem_all_tot=mem_all_tot+mem_all
      call myalloc_DYN() ; mem_all_tot=mem_all_tot+mem_all

#ifdef key_trc_dmp
!     needs Time_Manager
      call alloc_DTATRC(); mem_all_tot=mem_all_tot+mem_all
#endif
      call alloc_DIA()   ; mem_all_tot=mem_all_tot+mem_all

      call myalloc_BIO() ; mem_all_tot=mem_all_tot+mem_all
      call myalloc_SED() ; mem_all_tot=mem_all_tot+mem_all

      call myalloc_FN()  ; mem_all_tot=mem_all_tot+mem_all


      call MPI_ALLREDUCE(jpi, jpi_max, 1, MPI_INTEGER, MPI_MAX,MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(jpj, jpj_max, 1, MPI_INTEGER, MPI_MAX,MPI_COMM_WORLD, ierr)

      call myalloc_IO()  ; mem_all_tot=mem_all_tot+mem_all

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


      !call ave_init
      traIO     = trn ! initialization of Average  Matrix
      snIO      = 0.0
      tnIO      = 0.0
      vatmIO    = 0.0
      empIO     = 0.0
      qsrIO     = 0.0
      unIO      = 0.0
      bblxIO    = 0.0
      vnIO      = 0.0
      bblyIO    = 0.0
      wnIO      = 0.0
      avtIO     = 0.0
      tra_ppIO  = 0.0

      un        = 0.0
      vn        = 0.0
      wn        = 0.0
      call photo_init


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
        else
           if (lwp) write(*,*) 'Problems with start End. Program will stop. '
           STOP
      endif
! -----------------------------------------


        call getTimesteps(TimeStepStart, TimeStep__End)



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
