      SUBROUTINE forcings_PHYS(datestring)
!---------------------------------------------------------------------
!
!                       ROUTINE DTADYN
!                     ******************
!
!  PURPOSE :
!  ---------
!     Prepares dynamics and physics fields from an OPA8 simulation
!     for an off-line simulation for passive tracer
!                          =======
!
!   METHOD :
!   -------
!      calculates the position of DATA to read
!      READ DATA WHEN needed (example month changement)
!      computes slopes IF needed
!      interpolates DATA IF needed

!
!   EXTERNAL :                   prihre
!   --------
!
!----------------------------------------------------------------------
! parameters and commons
! ======================


       USE myalloc
       USE myalloc_mpp
       USE TIME_MANAGER
       IMPLICIT NONE

      character(LEN=17), INTENT(IN) ::  datestring

! local declarations
! ==================
      REAL(8) sec,zweigh
      integer Before, After
      INTEGER iper,iperm1,iswap


!     iper   : INTEGER for the period modulo ndtadyn
!     iperm1 : INTEGER for the precedent period  modulo ndtadyn
!     iswap  : indicator of swap of dynamic DATA array




       forcing_phys_partTime = MPI_WTIME()  ! cronometer-start



      sec=datestring2sec(DATEstring)
      call TimeInterpolation(sec,TC_FOR, BEFORE, AFTER, zweigh) ! 3.e-05 sec


      iperm1 = 1        ! li tengo buoni solo per la slopes
      iper = iperm1 + 1 !
      iswap  = 0


! ----------------------- INITIALISATION -------------
      IF (datestring.eq.DATESTART) then

          CALL LOAD_PHYS(TC_FOR%TimeStrings(TC_FOR%Before)) ! CALL dynrea(iperm1)
          call SLOPES(iperm1) ! Computes wdta (and slopes if key_trahdfiso)
          iswap = 1
          call swap_PHYS


        CALL LOAD_PHYS(TC_FOR%TimeStrings(TC_FOR%After)) !CALL dynrea(iper)
        call SLOPES(iper)! Computes wdta (and slopes if key_trahdfiso)

      ENDIF





! --------------------------------------------------------
! and now what we have to DO at every time step
! --------------------------------------------------------

! check the validity of the period in memory

      if (BEFORE.ne.TC_FOR%Before) then
         TC_FOR%Before = BEFORE
         TC_FOR%After  = AFTER

         call swap_PHYS
         iswap = 1


          CALL LOAD_PHYS(TC_FOR%TimeStrings(TC_FOR%After))
          CALL SLOPES(iper) !Computes wdta (and slopes if key_trahdfiso)

          IF(lwp) WRITE (numout,*) ' dynamics DATA READ for Time = ', TC_FOR%TimeStrings(TC_FOR%After)

!      ******* LOADED NEW FRAME *************
      END IF





! compute the DATA at the given time step

      SELECT CASE (nsptint)
           CASE (0)  !  ------- no time interpolation
!      we have to initialize DATA IF we have changed the period
              IF (iswap.eq.1) THEN
                 zweigh = 1.0
                 call ACTUALIZE_PHYS(zweigh)! initialize now fields with the NEW DATA READ
              END IF

          CASE (1) ! ------------linear interpolation ---------------

             call ACTUALIZE_PHYS(zweigh)



      END SELECT


! lb in any case, we need rhopn

!      jtask = 1; kt = 1 ; CALL eos(jtask,kt)




       forcing_phys_partTime = MPI_WTIME() - forcing_phys_partTime
       forcing_phys_TotTime  = forcing_phys_TotTime  + forcing_phys_partTime


      END SUBROUTINE forcings_PHYS

! ******************************************************
!     SUBROUTINE LOAD_PHYS(datestring)
!
!
! ******************************************************
       SUBROUTINE LOAD_PHYS(datestring)
! ======================
      USE calendar
      USE myalloc
      USE myalloc_mpp
      USE TIME_MANAGER

      IMPLICIT NONE

      CHARACTER(LEN=17), INTENT(IN) :: datestring
      LOGICAL B
      integer ji,jj,jk
! omp variables
            INTEGER :: mytid, ntids

#ifdef __OPENMP
            INTEGER ::  omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
            EXTERNAL :: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
#endif
      ! LOCAL
      character(LEN=30) nomefile


#ifdef __OPENMP
      ntids = omp_get_max_threads() ! take the number of threads
      mytid = -1000000
#else
      ntids = 1
      mytid = 0
#endif

      nomefile='FORCINGS/U19951206-12:00:00.nc'

! Starting I/O
! U  *********************************************************
      nomefile = 'FORCINGS/U'//datestring//'.nc'
      if(lwp) write(*,'(A,I4,A,A)') "LOAD_PHYS --> I am ", rank, " starting reading forcing fields from ", nomefile(1:30)
      call readnc_slice_float(nomefile,'vozocrtx',buf);

      DO jk=1,jpk,ntids!udta(:,:,:,2) = buf*umask;
!$omp parallel default(none) private(mytid,ji,jj) shared(jk,jpi,jpj,jpk,udta,umask,buf)
#ifdef __OPENMP
        mytid = omp_get_thread_num()  ! take the thread ID
#endif
         IF (jk+mytid <=jpk) then
            DO jj= 1,jpj
            DO ji= 1,jpi
              udta(ji,jj,jk+mytid,2) = buf(ji,jj,jk+mytid)*umask(ji,jj,jk+mytid)
            ENDDO
            ENDDO
         ENDIF
!$omp end parallel
      END DO



! V *********************************************************
      nomefile = 'FORCINGS/V'//datestring//'.nc'
      call readnc_slice_float(nomefile,'vomecrty',buf);

      DO jk=1,jpk,ntids!vdta(:,:,:,2) = buf*vmask;
!$omp parallel default(none) private(mytid,ji,jj) shared(jk,jpi,jpj,jpk,vdta,vmask,buf)
#ifdef __OPENMP
        mytid = omp_get_thread_num()  ! take the thread ID
#endif
         IF (jk+mytid <=jpk) then
            DO jj= 1,jpj
            DO ji= 1,jpi
              vdta(ji,jj,jk+mytid,2) = buf(ji,jj,jk+mytid)*vmask(ji,jj,jk+mytid)
            ENDDO
            ENDDO
         ENDIF
!$omp end parallel
      END DO


! W *********************************************************
      nomefile = 'FORCINGS/W'//datestring//'.nc'
      call readnc_slice_float(nomefile,'vovecrtz',buf);
      DO jk=1,jpk,ntids!wdta(:,:,:,2) = buf*tmask;
!$omp parallel default(none) private(mytid,ji,jj) shared(jk,jpi,jpj,jpk,wdta,tmask,buf)
#ifdef __OPENMP
        mytid = omp_get_thread_num()  ! take the thread ID
#endif
         IF (jk+mytid <=jpk) then
            DO jj= 1,jpj
            DO ji= 1,jpi
              wdta(ji,jj,jk+mytid,2) = buf(ji,jj,jk+mytid)*tmask(ji,jj,jk+mytid)
            ENDDO
            ENDDO
         ENDIF
!$omp end parallel
      END DO

      call readnc_slice_float(nomefile,'votkeavt',buf);
      DO jk=1,jpk,ntids!avtdta(:,:,:,2) = buf*tmask;
!$omp parallel default(none) private(mytid,ji,jj) shared(jk,jpi,jpj,jpk,avtdta,tmask,buf)
#ifdef __OPENMP
        mytid = omp_get_thread_num()  ! take the thread ID
#endif
         IF (jk+mytid <=jpk) then
            DO jj= 1,jpj
            DO ji= 1,jpi
              avtdta(ji,jj,jk+mytid,2) = buf(ji,jj,jk+mytid)*tmask(ji,jj,jk+mytid)
            ENDDO
            ENDDO
         ENDIF
!$omp end parallel
      END DO

! T *********************************************************
      nomefile = 'FORCINGS/T'//datestring//'.nc'
      call readnc_slice_float(nomefile,'votemper',buf);
      DO jk=1,jpk,ntids!tdta(:,:,:,2) = buf*tmask;
!$omp parallel default(none) private(mytid,ji,jj) shared(jk,jpi,jpj,jpk,tdta,tmask,buf)
#ifdef __OPENMP
        mytid = omp_get_thread_num()  ! take the thread ID
#endif
         IF (jk+mytid <=jpk) then
            DO jj= 1,jpj
            DO ji= 1,jpi
              tdta(ji,jj,jk+mytid,2) = buf(ji,jj,jk+mytid)*tmask(ji,jj,jk+mytid)
            ENDDO
            ENDDO
         ENDIF
!$omp end parallel
      END DO

      call readnc_slice_float(nomefile,'vosaline',buf);
      DO jk=1,jpk,ntids!sdta(:,:,:,2) = buf*tmask;
!$omp parallel default(none) private(mytid,ji,jj) shared(jk,jpi,jpj,jpk,sdta,tmask,buf)
#ifdef __OPENMP
        mytid = omp_get_thread_num()  ! take the thread ID
#endif
         IF (jk+mytid <=jpk) then
            DO jj= 1,jpj
            DO ji= 1,jpi
              sdta(ji,jj,jk+mytid,2) = buf(ji,jj,jk+mytid)*tmask(ji,jj,jk+mytid)
            ENDDO
            ENDDO
         ENDIF
!$omp end parallel
      END DO

      call readnc_slice_float_2d(nomefile,'sowindsp',buf2); flxdta(:,:,jpwind,2) = buf2*tmask(:,:,1);
      call readnc_slice_float_2d(nomefile,'soshfldo',buf2); flxdta(:,:,jpqsr ,2) = buf2*tmask(:,:,1);
                                                            flxdta(:,:,jpice ,2) = 0.
      call EXISTVAR(nomefile,'sowaflcd',B)
      if (B) then
      call readnc_slice_float_2d(nomefile,'sowaflcd',buf2); flxdta(:,:,jpemp ,2) = buf2*tmask(:,:,1);
      else
         if(lwp) write(*,*) 'Evaporation data not found. Forced to zero.'
         flxdta(:,:,jpemp ,2) = 0.
      endif


      bblxdta = 0.
      bblydta = 0.




#ifdef key_trahdfeiv
             udta(:,:,:,2) = udta(:,:,:,2) +zaeiu*umask
             vdta(:,:,:,2) = vdta(:,:,:,2) +zaeiv*vmask
             wdta(:,:,:,2) = wdta(:,:,:,2) +zaeiw*tmask
#endif



      END SUBROUTINE LOAD_PHYS





! ******************************************************
!     SUBROUTINE ACTUALIZE_PHYS(zweigh)
!     performs time interpolation
!     x(1)*(1-zweigh) + x(2)*zweigh
! ******************************************************
      SUBROUTINE ACTUALIZE_PHYS(zweigh)
         USE myalloc
         USE OPT_mem
         USE stafun
         IMPLICIT NONE
         REAL(8) zweigh, Umzweigh

         INTEGER ji,jj,jk
         INTEGER uk, uj      ! aux variables for OpenMP

      INTEGER :: mytid, ntids! omp variables

#ifdef __OPENMP
      INTEGER ::  omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
      EXTERNAL :: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
#endif


#ifdef __OPENMP
      ntids = omp_get_max_threads() ! take the number of threads
      mytid = -1000000

#else
      ntids = 1
      mytid = 0
#endif

      Umzweigh  = 1.0 - zweigh

          DO jk=1,jpk, ntids
!$omp parallel default(none) private(mytid,jj,ji)
!$omp&                       shared(jpi,jpj,jpk,jk,uk,ub,un,udta, vb,vn,vdta,wn,wdta,avt,avtdta,tn,tdta,sn,sdta,
!$omp&                              zweigh,Umzweigh,tmask)
#ifdef __OPENMP
         mytid = omp_get_thread_num()  ! take the thread ID
#endif
        if (mytid+jk.le.jpk) then

            uk = mytid+jk
            DO jj=1,jpj
              DO ji=1,jpi

                if (tmask(ji,jj,uk). NE . 0.0) then
                !ub(ji,jj,uk) = un(ji,jj,uk) ! serve ?
                !vb(ji,jj,uk) = vn(ji,jj,uk)
                un(ji,jj,uk) = (Umzweigh*  udta(ji,jj,uk,1) + zweigh*  udta(ji,jj,uk,2))
                vn(ji,jj,uk) = (Umzweigh*  vdta(ji,jj,uk,1) + zweigh*  vdta(ji,jj,uk,2))
                wn(ji,jj,uk) = (Umzweigh*  wdta(ji,jj,uk,1) + zweigh*  wdta(ji,jj,uk,2))
               avt(ji,jj,uk) = (Umzweigh*avtdta(ji,jj,uk,1) + zweigh*avtdta(ji,jj,uk,2))
#if defined key_off_tra
                tn(ji,jj,uk) = (Umzweigh*  tdta(ji,jj,uk,1) + zweigh*  tdta(ji,jj,uk,2))
                sn(ji,jj,uk) = (Umzweigh*  sdta(ji,jj,uk,1) + zweigh*  sdta(ji,jj,uk,2))
#endif
#if defined key_trahdfiso
              uslp(ji,jj,uk) = (Umzweigh *  uslpdta(ji,jj,uk,1)+ zweigh*  uslpdta(ji,jj,uk,2))
              vslp(ji,jj,uk) = (Umzweigh *  vslpdta(ji,jj,uk,1)+ zweigh*  vslpdta(ji,jj,uk,2))
             wslpi(ji,jj,uk) = (Umzweigh * wslpidta(ji,jj,uk,1)+ zweigh* wslpidta(ji,jj,uk,2))
             wslpj(ji,jj,uk) = (Umzweigh * wslpjdta(ji,jj,uk,1)+ zweigh* wslpjdta(ji,jj,uk,2))
#endif
#if  defined key_trahdfeiv 
              ugm(ji,jj,uk)  = (Umzweigh *   ugmdta(ji,jj,uk,1)+ zweigh*    ugmdta(ji,jj,uk,2))
              vgm(ji,jj,uk)  = (Umzweigh *   vgmdta(ji,jj,uk,1)+ zweigh*    vgmdta(ji,jj,uk,2))
              wgm(ji,jj,uk)  = (Umzweigh *   wgmdta(ji,jj,uk,1)+ zweigh*    wgmdta(ji,jj,uk,2))
#endif
              endif ! tmask
              END DO 
            END DO
      endif
!$omp end parallel
          END DO

          DO jj=1,jpj,ntids
!$omp parallel default(none) private(mytid,jk,ji)
!$omp&                       shared(jpi,jpj,jpk,jj,uj,flx,flxdta,bblx,bblxdta,bbly,bblydta,hmld,zmxldta,
!$omp&                              vatm,freeze,emp,qsr,jpwind,jpice,jpemp,jpqsr,zweigh, Umzweigh,jpflx)
#ifdef __OPENMP
         mytid = omp_get_thread_num()  ! take the thread ID
#endif
        uj = jj+mytid
        if (mytid+jj.le.jpj) then


            DO jk=1,jpflx
              DO ji=1,jpi
                flx(ji,uj,jk) = ( Umzweigh * flxdta(ji,uj,jk,1)+ zweigh     * flxdta(ji,uj,jk,2) )
              END DO 
            END DO


!          do jj=1,jpj
            do ji=1,jpi
              bblx(ji,uj)= Umzweigh * bblxdta(ji,uj,1)+zweigh    * bblxdta(ji,uj,2)
              bbly(ji,uj)= Umzweigh * bblydta(ji,uj,1)+zweigh    * bblydta(ji,uj,2)
              hmld(ji,uj)=(Umzweigh * zmxldta(ji,uj,1)+zweigh    * zmxldta(ji,uj,2))
            end do
!          end do

            DO ji=1,jpi
!                DO jj=1,jpj
                  vatm(ji,uj)   = flx(ji,uj,jpwind)
                  freeze(ji,uj) = flx(ji,uj,jpice)
                  emp(ji,uj)    = flx(ji,uj,jpemp)
                  qsr(ji,uj)    = flx(ji,uj,jpqsr)
            END DO

        endif
!$omp end parallel
       END DO


      END SUBROUTINE ACTUALIZE_PHYS



! *************************************************************
!      SUBROUTINE SWAP
! *    copia l'indice 2 nell'indice 1
! *************************************************************

      SUBROUTINE swap_PHYS
         USE myalloc
         IMPLICIT NONE
         INTEGER ji,jj,jk,jdepth
         INTEGER :: mytid, ntids! omp variables

#ifdef __OPENMP
      INTEGER ::  omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
      EXTERNAL :: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
#endif


#ifdef __OPENMP
      ntids = omp_get_max_threads() ! take the number of threads
      mytid = -1000000

#else
      ntids = 1
      mytid = 0
#endif

          DO jk=1,jpk,ntids
!$omp parallel default(None) private(mytid,ji,jj,jdepth) shared(jpi,jpj,jpk,jk,udta,vdta,wdta,avtdta,tdta,sdta)
#ifdef __OPENMP
         mytid = omp_get_thread_num()  ! take the thread ID
#endif
         jdepth=jk+mytid
         if (jdepth <= jpk) then

            DO jj=1,jpj
              DO ji=1,jpi
                  udta(ji,jj,jdepth,1) =  udta(ji,jj,jdepth,2)
                  vdta(ji,jj,jdepth,1) =  vdta(ji,jj,jdepth,2)
                  wdta(ji,jj,jdepth,1) =  wdta(ji,jj,jdepth,2)
                avtdta(ji,jj,jdepth,1) =avtdta(ji,jj,jdepth,2)
#if defined key_off_tra
                tdta(ji,jj,jdepth,1) = tdta(ji,jj,jdepth,2)
                sdta(ji,jj,jdepth,1) = sdta(ji,jj,jdepth,2)
#endif
#if defined key_trahdfiso
                 uslpdta(ji,jj,jdepth,1)=  uslpdta(ji,jj,jdepth,2)
                 vslpdta(ji,jj,jdepth,1)=  vslpdta(ji,jj,jdepth,2)
                wslpidta(ji,jj,jdepth,1)= wslpidta(ji,jj,jdepth,2)
                wslpjdta(ji,jj,jdepth,1)= wslpjdta(ji,jj,jdepth,2)
#endif
#if  defined key_trahdfeiv
                ugmdta(ji,jj,jdepth,1)=ugmdta(ji,jj,jdepth,2)
                vgmdta(ji,jj,jdepth,1)=vgmdta(ji,jj,jdepth,2)
                wgmdta(ji,jj,jdepth,1)=wgmdta(ji,jj,jdepth,2)
#endif
              END DO
            END DO
          ENDIF
!$omp end parallel

          END DO

          DO jk=1,jpflx,ntids
!$omp parallel default(None) private(mytid,ji,jj) shared(jk,jpi,jpj,jpflx,flxdta)
#ifdef __OPENMP
         mytid = omp_get_thread_num()  ! take the thread ID
#endif
           if (jk+mytid.le.jpflx) then
            DO ji=1,jpi
              DO jj=1,jpj
                flxdta(ji,jj,jk+mytid,1) = flxdta(ji,jj,jk+mytid,2)
              END DO
            END DO
           endif
!$omp end parallel
          END DO

          do jj=1,jpj
            do ji=1,jpi
              bblxdta(ji,jj,1) = bblxdta(ji,jj,2)
              bblydta(ji,jj,1) = bblydta(ji,jj,2)
              zmxldta(ji,jj,1) = zmxldta(ji,jj,2)
            end do
          end do


      END SUBROUTINE swap_PHYS

! *************************************************************
! SUBROUTINE SLOPES
! *
! *************************************************************

      SUBROUTINE SLOPES(iper)
         USE myalloc
         USE OPT_mem
         IMPLICIT NONE
         INTEGER iper
         INTEGER ji,jj,jk


#if defined key_off_tra
          DO jk=1,jpk
            DO jj=1,jpj
              DO ji=1,jpi
                tn(ji,jj,jk)=tdta(ji,jj,jk,2)
                sn(ji,jj,jk)=sdta(ji,jj,jk,2)
              END DO
            END DO
          END DO
#endif

!          CALL div(1,kt)
!          CALL wzv(1,kt)
#if defined key_trahdfiso
#if defined key_off_degrad
          CALL slprea(iper)
#else
          CALL eos()
          CALL bn2(   1,kt)
          CALL hdfslp(1,kt)
#endif
#endif
#if defined key_trahdfeiv && defined key_off_degrad
          CALL vgmrea(iper)
#endif

#if defined key_trahdfiso
          DO jk=1,jpk
            DO jj=1,jpj
              DO ji=1,jpi
                 uslpdta(ji,jj,jk,2)= uslp(ji,jj,jk)
                 vslpdta(ji,jj,jk,2)= vslp(ji,jj,jk)
                wslpidta(ji,jj,jk,2)=wslpi(ji,jj,jk)
                wslpjdta(ji,jj,jk,2)=wslpj(ji,jj,jk)
              END DO
            END DO
          END DO
#endif

      END SUBROUTINE SLOPES
