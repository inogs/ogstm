      SUBROUTINE step
!---------------------------------------------------------------------
!
!                       ROUTINE STEP
!                     ****************
!
!  PURPOSE :
!  ---------
!    Time loop of OPA
!
!   METHOD :
!   -------
!
!      READ the dynamics fiels and update for time step
!
!      Call passive tracer model trcstp
!                     and diagnostics trcdia
!                     and specific restart trcwri
!
!   INPUT :
!   -----
!
!
!   OUTPUT :            : no
!   ------
!
!   EXTERNAL :

!      trcstp, trcdia  passive tracers interface


       USE calendar
       USE myalloc
       USE myalloc_mpp
       USE TIME_MANAGER
       USE BC_mem
       IMPLICIT NONE


! local declarations
! ==================
      INTEGER TAU, indic
      REAL(8) stpparttime,stptottime
      character(LEN=17)  datestring, datemean, datefrom
      REAL(8) sec

!++++++++++++++++++++++++++++++c
!         Time  loop           c
!++++++++++++++++++++++++++++++c

        daytottime = 0.     !  Reset of the cronometer
        dtadyntottime = 0.
        dtatrctottime = 0.
        trcdiatottime = 0.
        trcstptottime = 0.
        trcwritottime = 0.
        trcdiatottime = 0.
        trcsmstottime = 0.
        trcopttottime = 0.
        trcbiotottime = 0.
        trcadvtottime = 0.
        trcdmptottime = 0.
        trclaphdftottime = 0.
        trcbilaphdftottime = 0.
        trcbbltottime = 0.
        trczdftottime = 0.
        trcnxttottime = 0.
        trcstptottime = 0.
        stptottime = 0.



      indic=1  ! 1 for indic to initialize output files (first call)

      if (IsStartBackup) then
         TauAVEfrom = datestringToTAU(BKPdatefrom)
      else
         TauAVEfrom = TimeStepStart -1
      endif



      DO TAU = TimeStepStart, TimeStep__End

         call tau2datestring(TAU, DATEstring)

         NOW_datestring = DATEstring ! update time manager module

         sec=datestring2sec(DATEstring)

         call yearly(DATEstring)

         stpparttime = MPI_WTIME()  ! stop cronomether
         if(lwp) write(*,*) "  step ----------------- Starting timestep = ",TAU,' time :', DATEstring


        if (IsaRestart(DATEstring)) then
            CALL trcwri(DATEstring) ! writes the restart files

         if(lwp) write(*,*) 'step backup = ', (.not.IsAnAveDump(DATEstring))

            if (.not.IsAnAveDump(DATEstring)) then ! backup conditions
               call tau2datestring(TauAVEfrom+1, datefrom)
               CALL trcdia(datestring, datefrom, datestring)
            endif

         endif









! For offline simulation READ DATA or precalculalted dynamics fields
! ------------------------------------------------------------------


      CALL forcings_PHYS(DATEstring)     ! dtadyn
      CALL forcings_KEXT(datestring)
      CALL bc_gib       (DATEstring)     ! CALL dtatrc(istp,0)! Gibraltar strait BC
      CALL bc_tin       (DATEstring)     ! CALL dtatrc(istp,1)
      CALL bc_atm       (DATEstring)     ! CALL dtatrc(istp,2)
      CALL eos          ()               ! Water density
! Call Passive tracer model between synchronization for small parallelisation
        CALL trcstp(TAU) ! se commento questo non fa calcoli
        call trcave
!        Diagnostics:
!        ------------

      if (IsAnAveDump(DATEstring)) then
         call MIDDLEDATE(TauAVEfrom, TAU, datemean)

         call tau2datestring(TauAVEfrom+1, datefrom)
         if (IsStartBackup) datefrom = BKPdatefrom ! overwrite

         CALL trcdia(datemean, datefrom, datestring)
         TauAVEfrom=TAU

         IsStartBackup = .false.
      endif






       stpparttime = MPI_WTIME() - stpparttime
       stptottime  = stptottime + stpparttime

! OPA TEMPORIZATION
       IF (TAU.GT.TimeStepStart) THEN
        IF( mod( TAU, nwritetrc ).EQ.0) THEN
           write(*,*) "************* OPA TEMPORIZATION ************"
           write(*,*) "              Iteration",TAU
           write(*,*) "routine******time_tot*********time_ave*****"
           write(*,*) "dtadyn____", dtadyntottime,"  ",dtadyntottime/(TAU - TimeStepStart +1), " rank->", rank
           write(*,*) "trcopt____", trcopttottime,"  ",trcopttottime/(TAU - TimeStepStart +1), " rank->", rank
           write(*,*) "trcsms____", trcsmstottime,"  ",trcsmstottime/(TAU - TimeStepStart +1), " rank->", rank
           write(*,*) "trcadv____", trcadvtottime,"  ",trcadvtottime/(TAU - TimeStepStart +1), " rank->", rank
           write(*,*) "trcbil____", trcbilaphdftottime,"  ",trcbilaphdftottime/(TAU - TimeStepStart +1), " rank->", rank
           write(*,*) "trcbbl____", trcbbltottime,"  ",trcbbltottime/(TAU - TimeStepStart +1), " rank->", rank
           write(*,*) "trczdf____", trczdftottime,"  ",trczdftottime/(TAU - TimeStepStart +1), " rank->", rank
           write(*,*) "trcnxt____", trcnxttottime,"  ",trcnxttottime/(TAU - TimeStepStart +1), " rank->", rank
           write(*,*) "trcstp____", trcstptottime,"  ",trcstptottime/(TAU - TimeStepStart +1), " rank->", rank
           write(*,*) "trcwri____", trcwritottime,"  ",trcwritottime/(TAU - TimeStepStart +1), " rank->", rank
           write(*,*) "trcdia____", trcdiatottime,"  ",trcdiatottime/(TAU - TimeStepStart +1), " rank->", rank
           write(*,*) "stp_______", stptottime   ,"  ",stptottime   /(TAU - TimeStepStart +1)
           write(*,*) "*******************************************"
       ENDIF
      ENDIF


!+++++++++++++++++++++++++++++c
!      End of time loop       c
!+++++++++++++++++++++++++++++c

      END DO  


      END SUBROUTINE step
