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
       USE IO_mem, only: ave_counter
       IMPLICIT NONE


! local declarations
! ==================
      INTEGER TAU, indic

      character(LEN=17)  datestring, datemean, datefrom
      REAL(8) sec
      LOGICAL B, isFIRST

!++++++++++++++++++++++++++++++c
!         Time  loop           c
!++++++++++++++++++++++++++++++c


      isFIRST=.true.

      indic=1  ! 1 for indic to initialize output files (first call)

      if (IsStartBackup) then
         TauAVEfrom = datestringToTAU(BKPdatefrom)
      else
         TauAVEfrom = TimeStepStart! -1
      endif



      DO TAU = TimeStepStart, TimeStep__End

         stpparttime = MPI_WTIME()  ! stop cronomether

         call tau2datestring(TAU, DATEstring)
         sec=datestring2sec(DATEstring)

         NOW_datestring = DATEstring ! update time manager module
         NOW_sec        = sec


         call yearly(DATEstring) ! Performs daily updates
         call daily(DATEstring)  ! Performs daily updates


         if(lwp) write(numout,'(A,I8,A,A)') "  step ----------------- Starting timestep = ",TAU,' time :', DATEstring
         if(lwp) write(*     ,'(A,I8,A,A)') "  step ----------------- Starting timestep = ",TAU,' time :', DATEstring


        if (IsaRestart(DATEstring)) then
            CALL trcwri(DATEstring) ! writes the restart files

            if(lwp) write(*,*) 'step backup = ', (.not.IsAnAveDump(DATEstring))

            if (.not.IsAnAveDump(DATEstring)) then ! backup conditions
               call tau2datestring(TauAVEfrom, datefrom)  ! il TauAVEfrom+1 mi causava la data dell'ave alle 12:15
               CALL trcdia(datestring, datefrom, datestring)
            endif
         B = writeTemporization("trcdia____", trcdiatottime)
         B = writeTemporization("trcwri____", trcwritottime)
         endif









! For offline simulation READ DATA or precalculalted dynamics fields
! ------------------------------------------------------------------


      CALL forcings_PHYS(DATEstring)     ! dtadyn
      CALL forcings_KEXT(datestring)
      CALL bc_gib       (DATEstring)     ! CALL dtatrc(istp,0)! Gibraltar strait BC
      CALL bc_tin       (DATEstring)     ! CALL dtatrc(istp,1)
      CALL bc_atm       (DATEstring)     ! CALL dtatrc(istp,2)
      CALL bc_co2       (DATEstring)
      CALL eos          ()               ! Water density
! Call Passive tracer model between synchronization for small parallelisation
        CALL trcstp(TAU) ! se commento questo non fa calcoli
        call trcave
        ave_counter = ave_counter +1  ! incrementing our counter
!        Diagnostics:
!        ------------

      if (IsAnAveDump(DATEstring)) then
         call MIDDLEDATE(TauAVEfrom, TAU, datemean)

         call tau2datestring(TauAVEfrom, datefrom)
         if (IsStartBackup) datefrom = BKPdatefrom ! overwrite

         CALL trcdia(datemean, datefrom, datestring)
         TauAVEfrom=TAU

         ave_counter = 0   !  reset the counter
         IsStartBackup = .false.
         B = writeTemporization("trcdia____", trcdiatottime)
      endif






       stpparttime = MPI_WTIME() - stpparttime
       stptottime  = stptottime  + stpparttime

! OPA TEMPORIZATION
       IF (TAU.GT.TimeStepStart) THEN
        IF( mod( TAU, nwritetrc ).EQ.0) THEN
           write(*,*) "************* OPA TEMPORIZATION ************"
           write(*,*) "              Iteration",TAU
           write(*,*) "routine******time_tot*********time_ave*****"
           B = writeTemporization("forPhys___", forcing_phys_TotTime)
           B = writeTemporization("forKext___", forcing_kext_TotTime)
           B = writeTemporization("bcCO2_____", bc_co2_TotTime)
           B = writeTemporization("bcTIN_____", bc_tin_TotTime)
           B = writeTemporization("bcATM_____", bc_atm_TotTime)
           B = writeTemporization("bcGIB_____", bc_gib_TotTime)
           B = writeTemporization("density___", density_TotTime)
           B = writeTemporization("averaging_", ave_TotTime   )
           B = writeTemporization("trcopt____", trcopttottime)
           B = writeTemporization("trcbio____", BIOtottime)
           B = writeTemporization("trcadv____", trcadvtottime)
           B = writeTemporization("trcdmp____", trcdmptottime)
           B = writeTemporization("trcbil____", trcbilaphdftottime)
           B = writeTemporization("trcbbl____", trcbbltottime)
           B = writeTemporization("trcsbc____", trcsbctottime)
           B = writeTemporization("trcsms____", trcsmstottime)
           B = writeTemporization("trczdf____", trczdftottime)
           B = writeTemporization("snutel____", snuteltottime)
           B = writeTemporization("check_____", checkVtottime)
           B = writeTemporization("trcnxt____", trcnxttottime)
           B = writeTemporization("trcstp____", trcstptottime)


           B = writeTemporization("flxdump___",flx_TotTime  )
           B = writeTemporization("stp_______", stptottime  )

           call reset_Timers()
           write(*,*) "*******************************************"
       ENDIF
      ENDIF


!+++++++++++++++++++++++++++++c
!      End of time loop       c
!+++++++++++++++++++++++++++++c

      END DO  

      CONTAINS

      LOGICAL FUNCTION writeTemporization(string, elapsedtime)
      IMPLICIT NONE
      CHARACTER(LEN=*) string
      REAL(8) elapsedtime

      if (isFIRST) then
         write(*,250) string,elapsedtime,elapsedtime/(TAU-TimeStepStart +1)," rank->", rank
         isFirst=.false.
      else
         write(*,250) string,elapsedtime,elapsedtime/nwritetrc," rank->", rank
      endif
      writeTemporization = .true.
250   FORMAT (A , ES11.4 ,ES20.7 ,A20 , I3 )
      END FUNCTION writeTemporization

      END SUBROUTINE step
