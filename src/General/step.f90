MODULE module_step

      USE calendar
      USE myalloc
      USE TIME_MANAGER
      USE BC_mem
      USE IO_mem, only: ave_counter_1, ave_counter_2
      USE mpi
      USE mod_atmbc
      USE mod_cbc
      USE mod_gibbc
      USE mod_tinbc
      USE ogstm_mpi_module


 implicit NONE

 contains
 
 SUBROUTINE step
!---------------------------------------------------------------------
!
!                       ROUTINE STEP
!                     ****************
!
!  PURPOSE :
!  ---------
!    Time loop of ogstm
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

       IMPLICIT NONE


! local declarations
! ==================
      INTEGER TAU, indic

      character(LEN=17)  datestring, datemean, datefrom_1, datefrom_2
      double precision sec
      LOGICAL B, isFIRST
      double precision :: timess,timesw
      INTEGER :: jk,jj,ji,jn
!++++++++++++++++++++++++++++++c
!         Time  loop           c
!++++++++++++++++++++++++++++++c


      isFIRST=.true.

      indic=1  ! 1 for indic to initialize output files (first call)

      TauAVEfrom_1 = TimeStepStart 
       if (IsStartBackup_1) TauAVEfrom_1 = datestringToTAU(BKPdatefrom_1)
      TauAVEfrom_2 = TimeStepStart 
       if (IsStartBackup_2) TauAVEfrom_2 = datestringToTAU(BKPdatefrom_2)

      DO TAU = TimeStepStart, TimeStep__End

         stpparttime = MPI_WTIME()  ! stop cronomether
         timess = MPI_WTIME()
         call tau2datestring(TAU, DATEstring)
         sec=datestring2sec(DATEstring)

         NOW_datestring = DATEstring ! update time manager module
         NOW_sec        = sec


         call yearly(DATEstring) ! Performs yearly updates
         call daily(DATEstring)  ! Performs daily  updates


         if(lwp) write(numout,'(A,I8,A,A)') "step ------------ Starting timestep = ",TAU,' time ',DATEstring
         if(lwp) write(*,'(A,I8,A,A)')      "step ------------ Starting timestep = ",TAU,' time ',DATEstring


        if (IsaRestart(DATEstring)) then
            CALL trcwri(DATEstring) ! writes the restart files



            if (.not.IsAnAveDump(DATEstring,1)) then ! backup conditions group 1
               call tau2datestring(TauAVEfrom_1, datefrom_1)
               CALL trcdia(datestring, datefrom_1, datestring,1)
            endif

            if (.not.IsAnAveDump(DATEstring,2)) then ! backup conditions group 2
               call tau2datestring(TauAVEfrom_2, datefrom_2)
               if (save_bkp_group2) CALL trcdia(datestring, datefrom_2, datestring,2)
            endif

         if (lwp) then
             B = writeTemporization("trcdia____", trcdiatottime)
             B = writeTemporization("trcwri____", trcwritottime)
         endif
         endif



! For offline simulation READ DATA or precalculalted dynamics fields
! ------------------------------------------------------------------


      CALL forcings_PHYS(DATEstring)
      CALL forcings_KEXT(datestring)
      CALL bc_gib       (DATEstring)     ! CALL dtatrc(istp,0)! Gibraltar strait BC
      CALL bc_tin       (DATEstring)     ! CALL dtatrc(istp,1)
      CALL bc_atm       (DATEstring)     ! CALL dtatrc(istp,2)
      CALL bc_co2       (DATEstring)
      CALL eos          ()               ! Water density

      timesw = MPI_WTIME()

      if (IsAnAveDump(DATEstring,1)) then
         call MIDDLEDATE(TauAVEfrom_1, TAU, datemean)

         call tau2datestring(TauAVEfrom_1, datefrom_1)
         if (IsStartBackup_1) datefrom_1 = BKPdatefrom_1 ! overwrite
         CALL trcdia(datemean, datefrom_1, datestring,1)
         TauAVEfrom_1    = TAU
         ave_counter_1   = 0   !  reset the counter
         IsStartBackup_1 = .false.
        if (lwp)  B = writeTemporization("trcdia____", trcdiatottime)
      endif

      if (IsAnAveDump(DATEstring,2)) then
         call MIDDLEDATE(TauAVEfrom_2, TAU, datemean)

         call tau2datestring(TauAVEfrom_2, datefrom_2)
         if (IsStartBackup_2) datefrom_2 = BKPdatefrom_2 ! overwrite
         CALL trcdia(datemean, datefrom_2, datestring,2)
         TauAVEfrom_2    = TAU
         ave_counter_2   = 0   !  reset the counter
         IsStartBackup_2 = .false.
         if (lwp) B = writeTemporization("trcdia____", trcdiatottime)
      endif

      timesw = MPI_WTIME() - timesw



! Call Passive tracer model between synchronization for small parallelisation
        CALL trcstp    ! se commento questo non fa calcoli

      ! OPEN(UNIT=10007, FILE='s7.txt', FORM='FORMATTED')
      ! DO jn=1,jptra; DO jk = 1,jpk; DO jj = 1,jpj; DO ji = 1,jpi;
      ! WRITE(10007,200),'S7',jn,jk,jj,ji,trn(jk,jj,ji,jn)
      ! ENDDO;ENDDO;ENDDO;ENDDO;CLOSE(10007)

        call trcave
        ave_counter_1 = ave_counter_1 +1  ! incrementing our counters
        ave_counter_2 = ave_counter_2 +1

!       OPEN(UNIT=10008, FILE='s8.txt', FORM='FORMATTED')
!       DO jn=1,jptra; DO jk = 1,jpk; DO jj = 1,jpj; DO ji = 1,jpi;
!       WRITE(10008,200),'S8',jn,jk,jj,ji,trn(jk,jj,ji,jn)
!       ENDDO;ENDDO;ENDDO;ENDDO;CLOSE(10008)
! 200     FORMAT(' ',A3,I4,I4,I4,I4,D30.23)
!        STOP
       stpparttime = MPI_WTIME() - stpparttime
       stptottime  = stptottime  + stpparttime
       timess = MPI_WTIME() - timess
       if (lwp) then
       print *,"TIME STEP s =",timess,stptottime
       endif
! OGSTM TEMPORIZATION
       IF (TAU.GT.TimeStepStart) THEN
        IF( mod( TAU, nwritetrc ).EQ.0) THEN
           if (lwp) then
               write(*,*) "************* OGSTM TEMPORIZATION *********"
               write(*,*) "              Iteration",TAU
               write(*,*) "routine******time_tot*********time_ave*****"
           endif
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
      double precision elapsedtime

      if (isFIRST) then
         write(*,250) string,elapsedtime,elapsedtime/(TAU-TimeStepStart +1)," myrank->", myrank
         isFirst=.false.
      else
         write(*,250) string,elapsedtime,elapsedtime/nwritetrc," myrank->", myrank
      endif
      writeTemporization = .true.
250   FORMAT (A , ES11.4 ,ES20.7 ,A20 , I3 )
      END FUNCTION writeTemporization

      END SUBROUTINE step


      SUBROUTINE trcstp
!---------------------------------------------------------------------
!
!                       ROUTINE trcstp
!                     *****************
!
!  PURPOSE :
!  ---------
!	time loop of ogstm for passive tracer
!
!   METHOD :
!   -------
!      compute the well/spring evolution
!
!      compute the time evolution of tracers concentration
!         with advection
!         with horizontal diffusion
!         with surface boundary condition
!         with IMPLICIT vertical diffusion

       ! epascolo USE myalloc_mpp
       IMPLICIT NONE
      integer jn,jk,ji,jj
      double precision :: timetrczdf
      trcstpparttime = MPI_WTIME() ! cronometer-start

      CALL trcadv ! tracers: advection

#    if defined key_trc_dmp
      CALL trcdmp ! tracers: damping for passive tracers
#    endif

! tracers: horizontal diffusion IF namelist flags are activated
! -----------------------------

      IF (lhdf)   CALL trchdf


! tracers: sink and source (must be  parallelized on vertical slab)
! ------------------------
      ! OPEN(UNIT=10001, FILE='s1.txt', FORM='FORMATTED')
      ! DO jn=1,jptra; DO jk = 1,jpk; DO jj = 1,jpj; DO ji = 1,jpi;
      ! WRITE(10001,200),'S1',jn,jk,jj,ji,tra(jk,jj,ji,jn)
      ! ENDDO;ENDDO;ENDDO;ENDDO; CLOSE(10001)

      CALL trcsms

      ! OPEN(UNIT=10002, FILE='s2.txt', FORM='FORMATTED')
      ! DO jn=1,jptra; DO jk = 1,jpk; DO jj = 1,jpj; DO ji = 1,jpi;
      ! WRITE(10002,200),'S2',jn,jk,jj,ji,tra(jk,jj,ji,jn)
      ! ENDDO;ENDDO;ENDDO;ENDDO;CLOSE(10002)
      timetrczdf = MPI_WTIME()
      CALL trczdf ! tracers: vertical diffusion
      timetrczdf = MPI_WTIME() - timetrczdf
      print *,"TIME ZDF",timetrczdf
      ! OPEN(UNIT=10003, FILE='s3.txt', FORM='FORMATTED')
      ! DO jn=1,jptra; DO jk = 1,jpk; DO jj = 1,jpj; DO ji = 1,jpi;
      ! WRITE(10003,200),'S3',jn,jk,jj,ji,tra(jk,jj,ji,jn)
      ! ENDDO;ENDDO;ENDDO;ENDDO;CLOSE(10003)

      CALL snutel
      
      ! OPEN(UNIT=10004, FILE='s4.txt', FORM='FORMATTED')
      ! DO jn=1,jptra; DO jk = 1,jpk; DO jj = 1,jpj; DO ji = 1,jpi;
      ! WRITE(10004,200),'S4',jn,jk,jj,ji,tra(jk,jj,ji,jn)
      ! ENDDO;ENDDO;ENDDO;ENDDO;CLOSE(10004)

      CALL checkValues

      ! OPEN(UNIT=10005, FILE='s5.txt', FORM='FORMATTED')
      ! DO jn=1,jptra; DO jk = 1,jpk; DO jj = 1,jpj; DO ji = 1,jpi;
      ! WRITE(10005,200),'S5',jn,jk,jj,ji,tra(jk,jj,ji,jn)
      ! ENDDO;ENDDO;ENDDO;ENDDO;CLOSE(10005)

      CALL trcnxt ! tracers: fields at next time step

!       OPEN(UNIT=10006, FILE='s6.txt', FORM='FORMATTED')
!       DO jn=1,jptra; DO jk = 1,jpk; DO jj = 1,jpj; DO ji = 1,jpi;
!       WRITE(10006,200),'S6',jn,jk,jj,ji,tra(jk,jj,ji,jn)
!       ENDDO;ENDDO;ENDDO;ENDDO;CLOSE(10006)

! 200     FORMAT(' ',A3,I4,I4,I4,I4,D30.23)

      

      trcstpparttime = MPI_WTIME() - trcstpparttime ! cronometer-stop
      trcstptottime = trcstptottime + trcstpparttime

      END SUBROUTINE trcstp

end module