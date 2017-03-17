      SUBROUTINE forcings_PHYS(datestring)
!---------------------------------------------------------------------
!
!                       ROUTINE DTADYN
!                     ******************
!
!  PURPOSE :
!  ---------
!     Prepares dynamics and physics fields
!     for an off-line simulation for passive tracer
!                          =======
!
!   METHOD :
!   -------
!      calculates the position of DATA to read
!      READ DATA WHEN needed (example month changement)
!      interpolates DATA IF needed
!
!----------------------------------------------------------------------
! parameters and commons
! ======================

       USE myalloc
       USE TIME_MANAGER
       use mpi
       IMPLICIT NONE

      character(LEN=17), INTENT(IN) ::  datestring

! local declarations
! ==================
      double precision :: sec,zweigh
      integer :: Before, After
      integer :: iswap



!     iswap  : indicator of swap of dynamic DATA array

       forcing_phys_partTime = MPI_WTIME()  ! cronometer-start

      sec=datestring2sec(DATEstring)
      call TimeInterpolation(sec,TC_FOR, BEFORE, AFTER, zweigh) ! 3.e-05 sec
 
      iswap  = 0

! ----------------------- INITIALISATION -------------
      IF (datestring.eq.DATESTART) then

          CALL LOAD_PHYS(TC_FOR%TimeStrings(TC_FOR%Before)) ! CALL dynrea(iperm1)


          iswap = 1
          call swap_PHYS


        CALL LOAD_PHYS(TC_FOR%TimeStrings(TC_FOR%After)) !CALL dynrea(iper)




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
      USE TIME_MANAGER

      IMPLICIT NONE

      character(LEN=17), INTENT(IN) :: datestring
      LOGICAL :: B
      integer  :: jk,jj,ji
      ! LOCAL
      character(LEN=30) nomefile

      nomefile='FORCINGS/U19951206-12:00:00.nc'

! Starting I/O
! U  *********************************************************
      nomefile = 'FORCINGS/U'//datestring//'.nc'
      if(lwp) write(*,'(A,I4,A,A)') "LOAD_PHYS --> I am ", myrank, " starting reading forcing fields from ", nomefile(1:30)
      call readnc_slice_float(nomefile,'vozocrtx',buf)
      

                  DO ji= 1,jpi
            DO jj= 1,jpj
      DO jk= 1,jpk
                    udta(jk,jj,ji,2) = buf(jk,jj,ji)*umask(jk,jj,ji)
                  ENDDO
            ENDDO
      ENDDO

      call readnc_slice_float(nomefile,'e3u',buf)
      

                  DO ji= 1,jpi
            DO jj= 1,jpj
      DO jk=1,jpk
                    e3udta(jk,jj,ji,2) = buf(jk,jj,ji)*umask(jk,jj,ji)
                  ENDDO
            ENDDO
      ENDDO



! V *********************************************************
      nomefile = 'FORCINGS/V'//datestring//'.nc'
      call readnc_slice_float(nomefile,'vomecrty',buf)
      

                  DO ji= 1,jpi
            DO jj= 1,jpj
      DO jk=1,jpk
                    vdta(jk,jj,ji,2) = buf(jk,jj,ji)*vmask(jk,jj,ji)
                  ENDDO
            ENDDO
      ENDDO

      call readnc_slice_float(nomefile,'e3v',buf)
      

                  DO ji= 1,jpi
            DO jj= 1,jpj
      DO jk=1,jpk
              e3vdta(jk,jj,ji,2) = buf(jk,jj,ji)*vmask(jk,jj,ji)
                  ENDDO
            ENDDO
      ENDDO



! W *********************************************************


      nomefile = 'FORCINGS/W'//datestring//'.nc'

      call readnc_slice_float(nomefile,'vovecrtz',buf)
      
                        DO ji= 1,jpi
            DO jj= 1,jpj
      DO jk=1,jpk
                   wdta(jk,jj,ji,2) = buf(jk,jj,ji)*tmask(jk,jj,ji)
                  ENDDO
            ENDDO
      ENDDO


      call readnc_slice_float(nomefile,'votkeavt',buf)
      
                  DO ji= 1,jpi
            DO jj= 1,jpj
      DO jk=1,jpk
              avtdta(jk,jj,ji,2) = buf(jk,jj,ji)*tmask(jk,jj,ji)
                  ENDDO
            ENDDO
      ENDDO


      nomefile = 'FORCINGS/W'//datestring//'.nc'
      call readnc_slice_float(nomefile,'e3w',buf)
      
                  DO ji= 1,jpi
            DO jj= 1,jpj
      DO jk=1,jpk
              e3wdta(jk,jj,ji,2) = buf(jk,jj,ji)*tmask(jk,jj,ji)
                   ENDDO
            ENDDO
      ENDDO


! T *********************************************************
      nomefile = 'FORCINGS/T'//datestring//'.nc'
      call readnc_slice_float(nomefile,'votemper',buf)
      
                  DO ji= 1,jpi
            DO jj= 1,jpj
      DO jk=1,jpk
                  tdta(jk,jj,ji,2) = buf(jk,jj,ji)*tmask(jk,jj,ji)
                   ENDDO
            ENDDO
      ENDDO

      call readnc_slice_float(nomefile,'vosaline',buf)
      
                 DO ji= 1,jpi
            DO jj= 1,jpj
      DO jk=1,jpk
                   sdta(jk,jj,ji,2) = buf(jk,jj,ji)*tmask(jk,jj,ji)
                   ENDDO
            ENDDO
      ENDDO

      call readnc_slice_float(nomefile,'e3t',buf)
      
                 DO ji= 1,jpi
            DO jj= 1,jpj
      DO jk=1,jpk
              e3tdta(jk,jj,ji,2) = buf(jk,jj,ji)*tmask(jk,jj,ji)
                   ENDDO
            ENDDO
      ENDDO



      call readnc_slice_float_2d(nomefile,'sowindsp',buf2)
      flxdta(:,:,jpwind,2) = buf2*tmask(1,:,:)
      call readnc_slice_float_2d(nomefile,'soshfldo',buf2)
      flxdta(:,:,jpqsr ,2) = buf2*tmask(1,:,:)
      flxdta(:,:,jpice ,2) = 0.
      call EXISTVAR(nomefile,'sowaflcd',B)
      if (B) then
            call readnc_slice_float_2d(nomefile,'sowaflcd',buf2)
            flxdta(:,:,jpemp ,2) = buf2*tmask(1,:,:)
      else
         if(lwp) write(*,*) 'Evaporation data not found. Forced to zero.'
         flxdta(:,:,jpemp ,2) = 0.
      endif

      call EXISTVAR(nomefile,'sossheiu',B)
      if (B) then
         call readnc_slice_float_2d(nomefile,'sossheiu',buf2)
           DO ji=1,jpi
         DO jj=1,jpj
            if (umask(1,jj,ji) .EQ. 1.) then  
                  flxdta(jj,ji,8 ,2) = buf2(jj,ji)
            endif      
           END DO
         END DO
      e3u(1,:,:) = flxdta(:,:,8 ,2)
      
      else
!     Do nothing leave the init value --> domrea
      endif
!!!!!!!!!!!!!!!!!
! epascolo warning
      call EXISTVAR(nomefile,'sossheiv',B)
      if (B) then
         call readnc_slice_float_2d(nomefile,'sossheiv',buf2)
            DO ji=1,jpi
         DO jj=1,jpj
               if (vmask(1,jj,ji) .EQ. 1.)  flxdta(jj,ji,9 ,2) = buf2(jj,ji)
      
            END DO
         END DO
      e3v(1,:,:) = flxdta(:,:,9 ,2)
      
      else
!     Do nothing leave the init value --> domrea
      endif
!!!!!!!!!!!!!!!!!

      call EXISTVAR(nomefile,'sossheit',B)
      if (B) then
         call readnc_slice_float_2d(nomefile,'sossheit',buf2)
            DO ji=1,jpi
         DO jj=1,jpj
               if (tmask(1,jj,ji) .EQ. 1.)  flxdta(jj,ji,10 ,2) = buf2(jj,ji)
            END DO
         END DO
      e3t(1,:,:) = flxdta(:,:,10 ,2)
      
      else
!     Do nothing leave the init value --> domrea
      endif
!!!!!!!!!!!!!!!!!


!     CALL div()               ! Horizontal divergence
!     CALL wzv()               ! vertical velocity

!        could be written for OpenMP
              DO ji=1,jpi
            DO jj=1,jpj
          DO jk=1,jpk
                tn(jk,jj,ji)=tdta(jk,jj,ji,2)
                sn(jk,jj,ji)=sdta(jk,jj,ji,2)
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
         double precision zweigh, Umzweigh

         INTEGER jk,jj,ji,jf
         INTEGER uk, uj      ! aux variables for OpenMP

   
      Umzweigh  = 1.0 - zweigh

!!!$omp parallel default(none) private(mytid,jj,ji,uk)
!!!$omp&                       shared(jpk,jpj,jpi,jk,ub,un,udta, vb,vn,vdta,wn,wdta,avt,avtdta,tn,tdta,sn,sdta,
!!!$omp&                       zweigh,Umzweigh,tmask,umask,vmask,e3u,e3udta,e3v,e3vdta,e3t,e3tdta,e3w,e3wdta,e3t_back)
          DO ji=1,jpi
          DO jj=1,jpj
          DO uk=1,jpk
                if (umask(uk,jj,ji) .eq. 1) then
                un(uk,jj,ji)  = (Umzweigh*  udta(uk,jj,ji,1) + zweigh*  udta(uk,jj,ji,2))
                e3u(uk,jj,ji) = (Umzweigh*  e3udta(uk,jj,ji,1) + zweigh*  e3udta(uk,jj,ji,2))
                endif
          ENDDO
          ENDDO
          ENDDO

          DO ji=1,jpi
          DO jj=1,jpj
          DO uk=1,jpk
                if (vmask(uk,jj,ji) .eq. 1) then
                vn(uk,jj,ji)  = (Umzweigh*  vdta(uk,jj,ji,1) + zweigh*  vdta(uk,jj,ji,2))
                e3v(uk,jj,ji) = (Umzweigh*  e3vdta(uk,jj,ji,1) + zweigh*  e3vdta(uk,jj,ji,2))
                endif
          ENDDO
          ENDDO
          ENDDO

          DO ji=1,jpi
          DO jj=1,jpj
          DO uk=1,jpk
                if (tmask(uk,jj,ji) .eq.1) then
                 wn(uk,jj,ji) = (Umzweigh*  wdta(uk,jj,ji,1) + zweigh*  wdta(uk,jj,ji,2))
                avt(uk,jj,ji) = (Umzweigh*avtdta(uk,jj,ji,1) + zweigh*avtdta(uk,jj,ji,2))
                e3w(uk,jj,ji) = (Umzweigh*  e3wdta(uk,jj,ji,1) + zweigh*  e3wdta(uk,jj,ji,2))
                endif
          ENDDO
          ENDDO
          ENDDO
          DO ji=1,jpi
          DO jj=1,jpj
          DO uk=1,jpk
                if (tmask(uk,jj,ji) .eq.1) then
                 tn(uk,jj,ji) = (Umzweigh*  tdta(uk,jj,ji,1) + zweigh*  tdta(uk,jj,ji,2))
                 sn(uk,jj,ji) = (Umzweigh*  sdta(uk,jj,ji,1) + zweigh*  sdta(uk,jj,ji,2))
                endif
          ENDDO
          ENDDO
          ENDDO

          DO ji=1,jpi
          DO jj=1,jpj
          DO uk=1,jpk
                if (tmask(uk,jj,ji) .eq.1) then          
                e3t_back(uk,jj,ji) = e3t(uk,jj,ji)
                e3t(uk,jj,ji) = (Umzweigh*  e3tdta(uk,jj,ji,1) + zweigh*  e3tdta(uk,jj,ji,2))
                endif ! tmask
          END DO
          END DO
          END DO




!!!$omp parallel default(none) private(mytid,jk,uj,ji)
!!!$omp&                       shared(jpk,jpj,jpi,jj,flx,flxdta,
!!!$omp&                              vatm,freeze,emp,qsr,jpwind,jpice,jpemp,jpqsr,zweigh, Umzweigh,jpflx)

                              DO jf=1,jpflx
                        DO ji=1,jpi
                  DO uj=1,jpj

                flx(uj,ji,jf) = ( Umzweigh * flxdta(uj,ji,jf,1)+ zweigh * flxdta(uj,ji,jf,2) )
                !if(jf==3) write(*,200),ji,uj,flx(uj,ji,jf)
                        END DO
                  END DO
            END DO
!            STOP
!            200 FORMAT(' ',I4,I4,D30.23)

                  DO ji=1,jpi
            DO uj=1,jpj
                  vatm(uj,ji)   = flx(uj,ji,jpwind)
                  freeze(uj,ji) = flx(uj,ji,jpice)
                  emp(uj,ji)    = flx(uj,ji,jpemp)
                  qsr(uj,ji)    = flx(uj,ji,jpqsr)
!                 e3u(uj,ji,1)  = flx(uj,ji,8)
!                 e3v(uj,ji,1)  = flx(uj,ji,9)
!                 e3t(uj,ji,1)  = flx(uj,ji,10)
            END DO
       END DO


      END SUBROUTINE ACTUALIZE_PHYS



! *************************************************************
!     SUBROUTINE SWAP
!     copies index 2 in index 1
! **************************************************************

      SUBROUTINE swap_PHYS
         USE myalloc
         IMPLICIT NONE
         INTEGER :: jk,jj,ji,jdepth,jf
         


              DO ji=1,jpi
            DO jj=1,jpj
          DO jdepth=1,jpk

                  udta(jdepth,jj,ji,1)   =  udta(jdepth,jj,ji,2)
                  e3udta(jdepth,jj,ji,1) =  e3udta(jdepth,jj,ji,2)
                  vdta(jdepth,jj,ji,1)   =  vdta(jdepth,jj,ji,2)
                  e3vdta(jdepth,jj,ji,1) =  e3vdta(jdepth,jj,ji,2)
                  wdta(jdepth,jj,ji,1)   =  wdta(jdepth,jj,ji,2)
                  e3wdta(jdepth,jj,ji,1) =  e3wdta(jdepth,jj,ji,2)
                  avtdta(jdepth,jj,ji,1) = avtdta(jdepth,jj,ji,2)
                    tdta(jdepth,jj,ji,1) = tdta(jdepth,jj,ji,2)
                    sdta(jdepth,jj,ji,1) = sdta(jdepth,jj,ji,2)
                  e3tdta(jdepth,jj,ji,1) =  e3tdta(jdepth,jj,ji,2)

              END DO
            END DO
          END DO

                        DO jf=1,jpflx
                   DO ji=1,jpi
              DO jj=1,jpj
                flxdta(jj,ji,jf,1) = flxdta(jj,ji,jf,2)
              END DO
            END DO
          END DO



      END SUBROUTINE swap_PHYS
