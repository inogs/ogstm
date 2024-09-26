MODULE module_step

      USE calendar
      USE myalloc
      USE TIME_MANAGER
      USE BC_mem
      USE IO_mem, only: elapsed_time_1, elapsed_time_2
      USE mpi
      USE mod_atmbc
      USE mod_cbc


! ----------------------------------------------------------------------
!  BEGIN BC_REFACTORING SECTION
!  ---------------------------------------------------------------------

      use bc_set_mod

! ----------------------------------------------------------------------
!  END BC_REFACTORING SECTION
!  ---------------------------------------------------------------------

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

       ! XXX: to be removed
       use DIA_mem, only: diaflx,flx_ridxt
       use myalloc, only: tra,trb,e1t,e3t_back,e2t,e3t,e3w,umask,vmask,tmask,avt,ahtt
       use BIO_mem, only: ogstm_sediPI,ogstm_PH,ogstm_co2
       USE OPT_mem, only: kef
       USE SED_mem
       USE ADV_mem

       use simple_timer
       IMPLICIT NONE


! local declarations
! ==================
      INTEGER TAU, indic

      character(LEN=17)  datestring, datemean, datefrom_1, datefrom_2
      character(LEN=17)  date_aveforDA
      LOGICAL B, isFIRST
      INTEGER :: jk,jj,ji,jn
!++++++++++++++++++++++++++++++c
!         Time  loop           c
!++++++++++++++++++++++++++++++c


      isFIRST=.true.

       datefrom_1 =  DATESTART
       if (IsStartBackup_1) datefrom_1 = BKPdatefrom_1
       datefrom_2 =  DATESTART
       if (IsStartBackup_2) datefrom_2 = BKPdatefrom_2
       datestring =  DATESTART
      TAU = 0
      call tstart("step_total")

      DO WHILE (.not.ISOVERTIME(datestring))

         call tstart("step")
         call tstart("step_1")
         stpparttime = MPI_WTIME()  ! stop cronomether
         COMMON_DATESTRING = DATEstring

         call yearly(DATEstring) ! Performs yearly updates
         call daily(DATEstring)  ! Performs daily  updates


         if(lwp) write(numout,'(A,I8,A,A)') "step ------------ Starting timestep = ",TAU,' time ',DATEstring
         if(lwp) write(*,'(A,I8,A,A)')      "step ------------ Starting timestep = ",TAU,' time ',DATEstring
         call tstop("step_1")

         call tstart("restart")
        if (IsaRestart(DATEstring)) then
            CALL trcwri(DATEstring) ! writes the restart files


            if (AVE_FREQ1%N .gt.0) then              !  void 1.aveTimes -> no backup
            if (.not.IsAnAveDump(DATEstring,1)) then ! backup conditions group 1
               CALL trcdia(datestring, datefrom_1, datestring,1)
            endif
            endif

            if (AVE_FREQ2%N .gt.0) then
            if (.not.IsAnAveDump(DATEstring,2)) then ! backup conditions group 2
               if (save_bkp_group2) CALL trcdia(datestring, datefrom_2, datestring,2)
            endif
            endif

         if (lwp) then
             B = writeTemporization("trcdia____", trcdiatottime)
             B = writeTemporization("trcwri____", trcwritottime)
         endif
         endif
         call tstop("restart")



! For offline simulation READ DATA or precalculalted dynamics fields
! ------------------------------------------------------------------

      call tstart("forcing_phys")
      CALL forcings_PHYS(DATEstring)
      !$acc update host(un,vn,wn,tn,sn,avt)
      !$acc update host(e3u,e3v,e3w) if(IS_FREE_SURFACE)
      !$acc update host(e3t_back,e3t) if(IS_FREE_SURFACE)
      !$acc update host(vatm,emp,qsr)
      call tstop("forcing_phys")
      call tstart("forcing_kext")
      CALL forcings_KEXT(datestring)
      !$acc update device(kef)
      call tstop("forcing_kext")

! ----------------------------------------------------------------------
!  BEGIN BC_REFACTORING SECTION
!  ---------------------------------------------------------------------

      call tstart("boundaries%update")
      call boundaries%update(datestring)
      call tstop("boundaries%update")

! ----------------------------------------------------------------------
!  END BC_REFACTORING SECTION
!  ---------------------------------------------------------------------

      call tstart("bc_atm")
      CALL bc_atm       (DATEstring)     ! CALL dtatrc(istp,2)
      call tstop("bc_atm")
      call tstart("bc_co2")
      CALL bc_co2       (DATEstring)
      call tstop("bc_co2")
      call tstart("eos")
      CALL eos          ()               ! Water density
      call tstop("eos")

      !$acc update device(tra,tmask)
      !$acc update device(zaa,zbb,zcc,inv_eu,inv_ev,inv_et,big_fact_zaa,big_fact_zbb,big_fact_zcc,zbtr_arr,e1t,e2t,e3t,e1u,e2u,e3u,e1v,e2v,e3v,e3w,trn,advmask,flx_ridxt,diaflx) if(ladv)
      !$acc update device(atm,e3t) if(latmosph)
      !$acc update device(umask,vmask,trb,ahtt,diaflx,flx_ridxt) if(lhdf)
      !$acc update device(e3t,rhopn,emp,trn) if (lsbc)
      !$acc update device(qsr,mbathy,bfmmask,trn,DAY_LENGTH,vatm,rho,e3t,gdept,ogstm_PH,ogstm_co2) if(lbfm)
#if  defined key_trc_sed
      !$acc update device(sed_idx,diaflx,e3t,ogstm_sedipi,mbathy) if(lbfm)
#endif
      !$acc update device(e1t,diaflx,e3t_back,e2t,trb,e3t,e3w) if (lzdf)
      !$acc update device(trn,umask,vmask,tmask,highfreq_table,e3t,tra_DIA,tra_DIA_2d,vatm,emp,qsr,highfreq_table_dia,highfreq_table_dia2d)

      call tstart("dump_ave_1")
      !$acc update host(traIO,traIO_HIGH,snIO,tnIO,wnIO,avtIO,e3tIO,unIO,vnIO,vatmIO,empIO,qsrIO,tra_DIA_IO,tra_DIA_2d_IO,tra_DIA_IO_HIGH,tra_DIA_2d_IO_HIGH)&
      !$acc& if(IsAnAveDump(DATEstring,1) .or. IsAnAveDump(DATEstring,2))
      if (IsAnAveDump(DATEstring,1)) then
         call MIDDLEDATE(datefrom_1, DATEstring, datemean)
         CALL trcdia(datemean, datefrom_1, datestring,1)

         datefrom_1      = DATEstring
         elapsed_time_1  = 0.0  !  reset the time counter
         IsStartBackup_1 = .false.

        if (lwp)  B = writeTemporization("trcdia____", trcdiatottime)
      endif
      call tstop("dump_ave_1")

      call tstart("dump_ave_2")
      if (IsAnAveDump(DATEstring,2)) then
         call MIDDLEDATE(datefrom_2, DATEstring, datemean)
         CALL trcdia(datemean, datefrom_2, datestring,2)
         datefrom_2      = DATEstring
         elapsed_time_2  = 0.0  !  reset the time counter
         IsStartBackup_2 = .false.
         if (lwp) B = writeTemporization("trcdia____", trcdiatottime)
      endif
      call tstop("dump_ave_2")


#ifdef ExecDA
      call tstart("data_assim")
      if (IsaDataAssimilation(DATEstring)) then
        CALL mainAssimilation(DATEstring, datefrom_1)
         if (lwp) B = writeTemporization("DATA_ASSIMILATION____", DAparttime)
      endif
      call tstop("data_assim")
#endif

! Call Passive tracer model between synchronization for small parallelisation
        call tstart("trcstp")
        CALL trcstp    ! se commento questo non fa calcoli
        call tstop("trcstp")
        call tstart("trcave")
        call trcave
        call tstop("trcave")

        elapsed_time_1 = elapsed_time_1 + rdt
        elapsed_time_2 = elapsed_time_2 + rdt


       stpparttime = MPI_WTIME() - stpparttime
       stptottime  = stptottime  + stpparttime


! OGSTM TEMPORIZATION
       call tstart("temporization")
       IF (TAU.GT.0) THEN
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
       call tstop("temporization")


!+++++++++++++++++++++++++++++c
!      End of time loop       c
!+++++++++++++++++++++++++++++c
      datestring = UPDATE_TIMESTRING(datestring, rdt)
      TAU = TAU + 1
         call tstop("step")

      !$acc update host(trb,trn,tra)
      !$acc update host(diaflx) if(lhdf)
      !$acc update host(zaa,zbb,zcc,inv_eu,inv_ev,inv_et,big_fact_zaa,big_fact_zbb,big_fact_zcc,zbtr_arr,diaflx) if(ladv)
      !$acc update host(tra_DIA,tra_DIA_2d,ogstm_sediPI,ogstm_PH) if(lbfm)
#if  defined key_trc_sed
      !$acc update host(diaflx,zwork) if(lbfm)
#endif
      !$acc update host(diaflx) if (lzdf)
      !$acc update host(tra_DIA_2d)
      !$acc update host(zwork) if(lbfm)
      !$acc update host(ogstm_sediPI) if(lbfm)
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

      use simple_timer
       IMPLICIT NONE
      integer jn,jk,ji,jj
      trcstpparttime = MPI_WTIME() ! cronometer-start

      call tstart("trcadv")

      IF (ladv) CALL trcadv ! tracers advection

      call tstop("trcadv")

#    if defined key_trc_dmp
      call tstart("trcdmp")
      CALL trcdmp ! tracers: damping for passive tracerstrcstp
      call tstop("trcdmp")

! ----------------------------------------------------------------------
!  BEGIN BC_REFACTORING SECTION
!  ---------------------------------------------------------------------

      call tstart("boundaries%apply")
      call boundaries%apply(e3t, trb, tra)
      call tstop("boundaries%apply")

! ----------------------------------------------------------------------
!  END BC_REFACTORING SECTION
!  ---------------------------------------------------------------------

#    endif

! tracers: horizontal diffusion IF namelist flags are activated
! -----------------------------

      call tstart("trchdf")
      IF (lhdf) CALL trchdf
      call tstop("trchdf")

! tracers: sink and source (must be  parallelized on vertical slab)
      call tstart("trcsbc")
      IF (lsbc) CALL trcsbc ! surface cell processes, default lsbc = False
      call tstop("trcsbc")

      call tstart("trcsms")

      IF (lbfm) CALL trcsms

      call tstop("trcsms")

      call tstart("trczdf")

      IF (lzdf) CALL trczdf ! tracers: vertical diffusion

      call tstop("trczdf")

      call tstart("snutel")
      IF (lsnu) CALL snutel
      call tstop("snutel")

      call boundaries%apply_dirichlet()

      ! CALL checkValues

      call tstart("trcnxt")
      CALL trcnxt ! tracers: fields at next time step
      call tstop("trcnxt")

      trcstpparttime = MPI_WTIME() - trcstpparttime ! cronometer-stop
      trcstptottime = trcstptottime + trcstpparttime

      END SUBROUTINE trcstp

end module
