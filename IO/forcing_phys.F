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
!      INTEGER kt   ! serve alla eos

!     INTEGER jtask
!     iper   : INTEGER for the period modulo ndtadyn
!     iperm1 : INTEGER for the precedent period  modulo ndtadyn
!     iswap  : indicator of swap of dynamic DATA array




       dtadynparttime = MPI_WTIME()  ! cronometer-start



      sec=datestring2sec(DATEstring)
      call TimeInterpolation(sec,TC_FOR, BEFORE, AFTER, zweigh)


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


       dtadynparttime = MPI_WTIME() - dtadynparttime
       dtadyntottime = dtadyntottime + dtadynparttime



      RETURN
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
      USE DYN_mem
      USE TIME_MANAGER

      IMPLICIT NONE

      CHARACTER(LEN=17), INTENT(IN) :: datestring


      ! LOCAL
      character(LEN=30) nomefile
      nomefile='FORCINGS/U19951206-12:00:00.nc'

! Starting I/O
! U  *********************************************************
      nomefile = 'FORCINGS/U'//datestring//'.nc'
      if(lwp) write(*,*) "LOAD_PHYS --> I am ", rank, " starting reading forcing fields from ", nomefile(1:30)
      call readnc_slice_float(nomefile,'vozocrtx',buf); udta(:,:,:,2) = buf*umask;

! V *********************************************************
      nomefile = 'FORCINGS/V'//datestring//'.nc'
      call readnc_slice_float(nomefile,'vomecrty',buf); vdta(:,:,:,2) = buf*vmask;


! W *********************************************************
      nomefile = 'FORCINGS/W'//datestring//'.nc'
      call readnc_slice_float(nomefile,'vovecrtz',buf);   wdta(:,:,:,2) = buf*tmask;
      call readnc_slice_float(nomefile,'votkeavt',buf); avtdta(:,:,:,2) = buf*tmask;

! T *********************************************************
      nomefile = 'FORCINGS/T'//datestring//'.nc'
      call readnc_slice_float(nomefile,'votemper',buf);   tdta(:,:,:,2) = buf*tmask;
      call readnc_slice_float(nomefile,'vosaline',buf);   sdta(:,:,:,2) = buf*tmask;

      call readnc_slice_float_2d(nomefile,'sowindsp',buf2); flxdta(:,:,jpwind,2) = buf2*tmask(:,:,1);
      call readnc_slice_float_2d(nomefile,'soshfldo',buf2); flxdta(:,:,jpqsr ,2) = buf2*tmask(:,:,1);
                                                            flxdta(:,:,jpice ,2) = 0.
                                                            flxdta(:,:,jpemp ,2) = 0.


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
         REAL(8) zweigh

         INTEGER ji,jj,jk
          DO jk=1,jpk
            DO jj=1,jpj
              DO ji=1,jpi
                ub(ji,jj,jk) = un(ji,jj,jk) ! serve ?
                vb(ji,jj,jk) = vn(ji,jj,jk)
                un(ji,jj,jk) = ((1. - zweigh)*  udta(ji,jj,jk,1) + zweigh*  udta(ji,jj,jk,2))
                vn(ji,jj,jk) = ((1. - zweigh)*  vdta(ji,jj,jk,1) + zweigh*  vdta(ji,jj,jk,2))
                wn(ji,jj,jk) = ((1. - zweigh)*  wdta(ji,jj,jk,1) + zweigh*  wdta(ji,jj,jk,2))
               avt(ji,jj,jk) = ((1. - zweigh)*avtdta(ji,jj,jk,1) + zweigh*avtdta(ji,jj,jk,2))
#if defined key_off_tra
                tn(ji,jj,jk) = ((1. - zweigh)*  tdta(ji,jj,jk,1) + zweigh*  tdta(ji,jj,jk,2))
                sn(ji,jj,jk) = ((1. - zweigh)*  sdta(ji,jj,jk,1) + zweigh*  sdta(ji,jj,jk,2))
#endif
#if defined key_trahdfiso
              uslp(ji,jj,jk) = ((1. - zweigh) *  uslpdta(ji,jj,jk,1)+ zweigh*  uslpdta(ji,jj,jk,2))
              vslp(ji,jj,jk) = ((1. - zweigh) *  vslpdta(ji,jj,jk,1)+ zweigh*  vslpdta(ji,jj,jk,2))
             wslpi(ji,jj,jk) = ((1. - zweigh) * wslpidta(ji,jj,jk,1)+ zweigh* wslpidta(ji,jj,jk,2))
             wslpj(ji,jj,jk) = ((1. - zweigh) * wslpjdta(ji,jj,jk,1)+ zweigh* wslpjdta(ji,jj,jk,2))
#endif
#if  defined key_trahdfeiv 
              ugm(ji,jj,jk)  = ((1. - zweigh) *   ugmdta(ji,jj,jk,1)+ zweigh*    ugmdta(ji,jj,jk,2))
              vgm(ji,jj,jk)  = ((1. - zweigh) *   vgmdta(ji,jj,jk,1)+ zweigh*    vgmdta(ji,jj,jk,2))
              wgm(ji,jj,jk)  = ((1. - zweigh) *   wgmdta(ji,jj,jk,1)+ zweigh*    wgmdta(ji,jj,jk,2))
#endif
              END DO 
            END DO
          END DO

          DO jk=1,jpflx
            DO ji=1,jpi
              DO jj=1,jpj
                flx(ji,jj,jk) = ( (1. - zweigh) * flxdta(ji,jj,jk,1)+ zweigh     * flxdta(ji,jj,jk,2) )
              END DO 
            END DO
          END DO

          do jj=1,jpj
            do ji=1,jpi
              bblx(ji,jj)= (1.-zweigh) * bblxdta(ji,jj,1)+zweigh    * bblxdta(ji,jj,2)
              bbly(ji,jj)= (1.-zweigh) * bblydta(ji,jj,1)+zweigh    * bblydta(ji,jj,2)
              hmld(ji,jj)=((1.-zweigh) * zmxldta(ji,jj,1)+zweigh    * zmxldta(ji,jj,2))
            end do
          end do

            DO ji=1,jpi
                DO jj=1,jpj
                  vatm(ji,jj)   = flx(ji,jj,jpwind)
                  freeze(ji,jj) = flx(ji,jj,jpice)
                  emp(ji,jj)    = flx(ji,jj,jpemp)
                  qsr(ji,jj)    = flx(ji,jj,jpqsr)
                END DO
              END DO


      END SUBROUTINE ACTUALIZE_PHYS



! *************************************************************
!      SUBROUTINE SWAP
! *    copia l'indice 2 nell'indice 1
! *************************************************************

      SUBROUTINE swap_PHYS
         USE myalloc
         IMPLICIT NONE
         INTEGER ji,jj,jk

          DO jk=1,jpk
            DO jj=1,jpj
              DO ji=1,jpi
                  udta(ji,jj,jk,1) =  udta(ji,jj,jk,2)
                  vdta(ji,jj,jk,1) =  vdta(ji,jj,jk,2)
                  wdta(ji,jj,jk,1) =  wdta(ji,jj,jk,2)
                avtdta(ji,jj,jk,1) =avtdta(ji,jj,jk,2)
#if defined key_off_tra
                tdta(ji,jj,jk,1) = tdta(ji,jj,jk,2)
                sdta(ji,jj,jk,1) = sdta(ji,jj,jk,2)
#endif
#if defined key_trahdfiso
                 uslpdta(ji,jj,jk,1)=  uslpdta(ji,jj,jk,2)
                 vslpdta(ji,jj,jk,1)=  vslpdta(ji,jj,jk,2)
                wslpidta(ji,jj,jk,1)= wslpidta(ji,jj,jk,2)
                wslpjdta(ji,jj,jk,1)= wslpjdta(ji,jj,jk,2)
#endif
#if  defined key_trahdfeiv
                ugmdta(ji,jj,jk,1)=ugmdta(ji,jj,jk,2)
                vgmdta(ji,jj,jk,1)=vgmdta(ji,jj,jk,2)
                wgmdta(ji,jj,jk,1)=wgmdta(ji,jj,jk,2)
#endif
              END DO
            END DO
          END DO

          DO jk=1,jpflx
            DO ji=1,jpi
              DO jj=1,jpj
                flxdta(ji,jj,jk,1) = flxdta(ji,jj,jk,2)
              END DO
            END DO
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
#    if defined key_off_degrad
          CALL slprea(iper)
#    else
          CALL eos()
          CALL bn2(   1,kt)
          CALL hdfslp(1,kt)
#    endif
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
