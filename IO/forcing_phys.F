
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
      INTEGER jtask
      INTEGER kt   ! serve alla eos

!     iper   : INTEGER for the period modulo ndtadyn
!     iperm1 : INTEGER for the precedent period  modulo ndtadyn
!     iswap  : indicator of swap of dynamic DATA array

! statement functions
! ===================

#include "stafun.h"


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
#include "stafun.h"

      ! LOCAL
      character*9  nomevar
      character*31 nomefile
      INTEGER jinew, jjnew, jista, jiend, jjsta, jjend
      INTEGER ji, jj, jk, ik
      REAL(8) zavmx
      LOGICAL condition

      nomevar= '123456789'
      nomefile='FORCINGS/U19951206-12:00:00.nc'

      nomevar(9:9)    = CHAR(0) ! once and for all
      nomefile(31:31) = CHAR(0)


! Starting I/O
! U  *********************************************************
      nomefile(1:30) = 'FORCINGS/U'//datestring//'.nc'
      if(lwp) write(*,*) "LOAD_PHYS --> I am ", rank, " starting reading forcing fields from ", nomefile(1:30)

      nomevar(1:8)='vozocrtx'; CALL ioogsnc_3du(nomefile,nomevar,1,jpiglo,jpjglo,jpk, zu_f)
!     nomevar(1:8)='sobblcox'; CALL ioogsnc_2du(nomefile,nomevar,1,jpiglo,jpjglo,  zbblx_f)
      zbblx_f=0.

! V *********************************************************
      nomefile(1:30) = 'FORCINGS/V'//datestring//'.nc'
      nomevar(1:8)='vomecrty'; CALL ioogsnc_3dv(nomefile,nomevar,1,jpiglo,jpjglo,jpk, zv_f)
!     nomevar(1:8)='sobblcoy'; CALL ioogsnc_2dv(nomefile,nomevar,1,jpiglo,jpjglo,zbbly)
      zbbly_f=0.

! W *********************************************************
      nomefile(1:30) = 'FORCINGS/W'//datestring//'.nc'
      nomevar(1:8)='vovecrtz';CALL ioogsnc_3dw(nomefile,nomevar,1,jpiglo,jpjglo,jpk,   zw_f)
      nomevar(1:8)='votkeavt';CALL ioogsnc_3dw(nomefile,nomevar,1,jpiglo,jpjglo,jpk, zavt_f)

! T *********************************************************
      nomefile(1:30) = 'FORCINGS/T'//datestring//'.nc'
      nomevar(1:8)='votemper'; CALL ioogsnc_3dt(nomefile,nomevar,1,jpiglo,jpjglo,jpk, zt_f)
      nomevar(1:8)='vosaline'; CALL ioogsnc_3dt(nomefile,nomevar,1,jpiglo,jpjglo,jpk, zs_f)
      nomevar(1:8)='sowindsp'; CALL ioogsnc_2dt(nomefile,nomevar,1,jpiglo,jpjglo,  zwind_f)
!     nomevar(1:8)='sowaflup'; CALL ioogsnc_2dt(nomefile,nomevar,1,jpiglo,jpjglo,   zemp_f)
       zemp_f = 0.
       zice_f = 0.
      nomevar(1:8) = 'soshfldo'; CALL ioogsnc_2dt(nomefile,nomevar,1,jpiglo,jpjglo, zqsr_f)



! In this way only one file has to be read, also with more than a cpu
      jista = nimpp
      jiend = MIN(jpi+nimpp-1, jpiglo)
      jjsta = njmpp
      jjend = MIN(jpj+njmpp-1, jpjglo)

! I change all the indexes of the loops
!      DO jk=1,jpk
!        DO jj=1,jpj
!          DO ji=1,jpi
      DO jk=1,jpk
        DO jj=jjsta, jjend
          DO ji=jista, jiend
            jinew = ji - nimpp + 1
            jjnew = jj - njmpp + 1
            udta(jinew,jjnew,jk,2)=REAL(zu_f(ji,jj,jk),8)*umask(jinew,jjnew,jk)
            vdta(jinew,jjnew,jk,2)=REAL(zv_f(ji,jj,jk),8)*vmask(jinew,jjnew,jk)
            wdta(jinew,jjnew,jk,2)=REAL(zw_f(ji,jj,jk),8)*tmask(jinew,jjnew,jk)

#ifdef key_trahdfeiv
             udta(jinew,jjnew,jk,2) = udta(jinew,jjnew,jk,2)+zaeiu(ji,jj,jk)*umask(jinew,jjnew,jk)
             vdta(jinew,jjnew,jk,2) = vdta(jinew,jjnew,jk,2)+zaeiv(ji,jj,jk)*vmask(jinew,jjnew,jk)
             wdta(jinew,jjnew,jk,2) = wdta(jinew,jjnew,jk,2)+zaeiw(ji,jj,jk)*tmask(jinew,jjnew,jk)
#endif

              tdta(jinew,jjnew,jk,2)=REAL(zt_f(  ji,jj,jk),8)*tmask(jinew,jjnew,jk)
              sdta(jinew,jjnew,jk,2)=REAL(zs_f(  ji,jj,jk),8)*tmask(jinew,jjnew,jk)
            avtdta(jinew,jjnew,jk,2)=REAL(zavt_f(ji,jj,jk),8)*tmask(jinew,jjnew,jk)
          END DO
        END DO
      END DO



      DO jj=jjsta, jjend
         DO ji=jista, jiend
            jinew = ji - nimpp + 1
            jjnew = jj - njmpp + 1
            flxdta(jinew,jjnew,jpwind,2)=REAL(zwind_f(ji,jj),8)*tmask(jinew,jjnew,1)
            flxdta(jinew,jjnew,jpice ,2)=min(1.D0,REAL(zice_f(ji,jj),8))*tmask(jinew,jjnew,1)
            flxdta(jinew,jjnew,jpemp, 2)=REAL(zemp_f(ji,jj),8)*tmask(jinew,jjnew,1)
            flxdta(jinew,jjnew,jpqsr, 2)=REAL(zqsr_f(ji,jj),8)*tmask(jinew,jjnew,1)
            ik = max( 31, 1 )
            bblxdta(jinew,jjnew,2)=max(0.D0,REAL(zbblx_f(ji,jj),8))*umask(jinew,jjnew,ik)
            bblydta(jinew,jjnew,2)=max(0.D0,REAL(zbbly_f(ji,jj),8))*vmask(jinew,jjnew,ik)
            if (bblxdta(jinew,jjnew,2).gt.2.) bblxdta(jinew,jjnew,2)=0.
            if (bblydta(jinew,jjnew,2).gt.2.) bblydta(jinew,jjnew,2)=0.

         END DO
      END DO


!  ---------------------------------------------------------
!  calcul de la profondeur de la couche de melange
!  en metre utile pour la formulation de la biologie de hamocc3
      zavmx = 5.0e-4
      imel=0
      DO jk=2,jpkm1
           DO jj=jjsta, jjend
             DO ji=jista, jiend
              jinew = ji - nimpp + 1
              jjnew = jj - njmpp + 1
              condition = zavt_f ( ji, jj, jk) .LE. zavmx.AND. imel(ji,jj) .EQ. 0
            IF (condition) THEN
                zmxldta (jinew,jjnew,2) = gdepw (jk+1)*tmask (jinew,jjnew,1)
                imel(ji,jj) = 1
            ENDIF
          END DO
        END DO
      END DO




      END SUBROUTINE LOAD_PHYS





! ******************************************************
!     SUBROUTINE ACTUALIZE_PHYS(zweigh)
!     performs time interpolation
!     x(1)*(1-zweigh) + x(2)*zweigh
! ******************************************************
      SUBROUTINE ACTUALIZE_PHYS(zweigh)
         USE myalloc
         USE OPT_mem
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
         INTEGER jtask


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
          CALL eos(   1,kt)
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
