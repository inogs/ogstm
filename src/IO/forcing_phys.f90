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
       ! epascolo USE myalloc_mpp
       USE TIME_MANAGER
       use mpi
       IMPLICIT NONE

      character(LEN=17), INTENT(IN) ::  datestring

! local declarations
! ==================
      REAL(8) sec,zweigh
      integer Before, After
      INTEGER iswap



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
      ! epascolo USE myalloc_mpp
      USE TIME_MANAGER

      IMPLICIT NONE

      CHARACTER(LEN=17), INTENT(IN) :: datestring
      LOGICAL B
      integer jk,jj,ji
! omp variables
            INTEGER :: mytid, ntids

#ifdef __OPENMP1
            INTEGER ::  omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
            EXTERNAL :: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
#endif
      ! LOCAL
      character(LEN=30) nomefile


#ifdef __OPENMP1
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
      if(lwp) write(*,'(A,I4,A,A)') "LOAD_PHYS --> I am ", myrank, " starting reading forcing fields from ", nomefile(1:30)
      call readnc_slice_float(nomefile,'vozocrtx',buf)
      

      DO jk=1,jpk,ntids!udta(:,:,:,2) = buf*umask
      
!!!$omp parallel default(none) private(mytid,jj,ji) shared(jk,jpk,jpj,jpi,udta,umask,buf)
#ifdef __OPENMP1
        mytid = omp_get_thread_num()  ! take the thread ID
#endif
         IF (jk+mytid <=jpk) then
            DO jj= 1,jpj
            DO ji= 1,jpi
              udta(jk,jj,ji+mytid,2) = buf(jk,jj,ji+mytid)*umask(jk,jj,ji+mytid)
            ENDDO
            ENDDO
         ENDIF
!!!$omp end parallel
      END DO

      call readnc_slice_float(nomefile,'e3u',buf)
      

      DO jk=1,jpk,ntids
!!!$omp parallel default(none) private(mytid,jj,ji) shared(jk,jpk,jpj,jpi,e3udta,umask,buf)
#ifdef __OPENMP1
        mytid = omp_get_thread_num()  ! take the thread ID
#endif
         IF (jk+mytid <=jpk) then
            DO jj= 1,jpj
            DO ji= 1,jpi
              e3udta(jk,jj,ji+mytid,2) = buf(jk,jj,ji+mytid)*umask(jk,jj,ji+mytid)
            ENDDO
            ENDDO
         ENDIF
!!!$omp end parallel
      END DO




! V *********************************************************
      nomefile = 'FORCINGS/V'//datestring//'.nc'
      call readnc_slice_float(nomefile,'vomecrty',buf)
      

      DO jk=1,jpk,ntids!vdta(:,:,:,2) = buf*vmask
      
!!!$omp parallel default(none) private(mytid,jj,ji) shared(jk,jpk,jpj,jpi,vdta,vmask,buf)
#ifdef __OPENMP1
        mytid = omp_get_thread_num()  ! take the thread ID
#endif
         IF (jk+mytid <=jpk) then
            DO jj= 1,jpj
            DO ji= 1,jpi
              vdta(jk,jj,ji+mytid,2) = buf(jk,jj,ji+mytid)*vmask(jk,jj,ji+mytid)
            ENDDO
            ENDDO
         ENDIF
!!!$omp end parallel
      END DO

      call readnc_slice_float(nomefile,'e3v',buf)
      

      DO jk=1,jpk,ntids
!!!$omp parallel default(none) private(mytid,jj,ji) shared(jk,jpk,jpj,jpi,e3vdta,vmask,buf)
#ifdef __OPENMP1
        mytid = omp_get_thread_num()  ! take the thread ID
#endif
         IF (jk+mytid <=jpk) then
            DO jj= 1,jpj
            DO ji= 1,jpi
              e3vdta(jk,jj,ji+mytid,2) = buf(jk,jj,ji+mytid)*vmask(jk,jj,ji+mytid)
            ENDDO
            ENDDO
         ENDIF
!!!$omp end parallel
      END DO



! W *********************************************************


      nomefile = 'FORCINGS/W'//datestring//'.nc'

      call readnc_slice_float(nomefile,'vovecrtz',buf)
      
      DO jk=1,jpk,ntids!wdta(:,:,:,2) = buf*tmask
      
!!!$omp parallel default(none) private(mytid,jj,ji) shared(jk,jpk,jpj,jpi,avtdta,tmask,buf,wdta)
#ifdef __OPENMP1
        mytid = omp_get_thread_num()  ! take the thread ID
#endif
         IF (jk+mytid <=jpk) then
            DO jj= 1,jpj
            DO ji= 1,jpi
              wdta(jk,jj,ji+mytid,2) = buf(jk,jj,ji+mytid)*tmask(jk,jj,ji+mytid)
            ENDDO
            ENDDO
         ENDIF
!!!$omp end parallel
      END DO


      call readnc_slice_float(nomefile,'votkeavt',buf)
      
      DO jk=1,jpk,ntids!avtdta(:,:,:,2) = buf*tmask
      
!!!$omp parallel default(none) private(mytid,jj,ji) shared(jk,jpk,jpj,jpi,avtdta,tmask,buf)
#ifdef __OPENMP1
        mytid = omp_get_thread_num()  ! take the thread ID
#endif
         IF (jk+mytid <=jpk) then
            DO jj= 1,jpj
            DO ji= 1,jpi
              avtdta(jk,jj,ji+mytid,2) = buf(jk,jj,ji+mytid)*tmask(jk,jj,ji+mytid)
            ENDDO
            ENDDO
         ENDIF
!!!$omp end parallel
      END DO


      nomefile = 'FORCINGS/W'//datestring//'.nc'
      call readnc_slice_float(nomefile,'e3w',buf)
      
      DO jk=1,jpk,ntids
!!!$omp parallel default(none) private(mytid,jj,ji) shared(jk,jpk,jpj,jpi,e3wdta,tmask,buf)
#ifdef __OPENMP1
        mytid = omp_get_thread_num()  ! take the thread ID
#endif
         IF (jk+mytid <=jpk) then
            DO jj= 1,jpj
            DO ji= 1,jpi
              e3wdta(jk,jj,ji+mytid,2) = buf(jk,jj,ji+mytid)*tmask(jk,jj,ji+mytid)
            ENDDO
            ENDDO
         ENDIF
!!!$omp end parallel
      END DO


! T *********************************************************
      nomefile = 'FORCINGS/T'//datestring//'.nc'
      call readnc_slice_float(nomefile,'votemper',buf)
      
      DO jk=1,jpk,ntids!tdta(:,:,:,2) = buf*tmask
      
!!!$omp parallel default(none) private(mytid,jj,ji) shared(jk,jpk,jpj,jpi,tdta,tmask,buf)
#ifdef __OPENMP1
        mytid = omp_get_thread_num()  ! take the thread ID
#endif
         IF (jk+mytid <=jpk) then
            DO jj= 1,jpj
            DO ji= 1,jpi
              tdta(jk,jj,ji+mytid,2) = buf(jk,jj,ji+mytid)*tmask(jk,jj,ji+mytid)
            ENDDO
            ENDDO
         ENDIF
!!!$omp end parallel
      END DO

      call readnc_slice_float(nomefile,'vosaline',buf)
      
      DO jk=1,jpk,ntids!sdta(:,:,:,2) = buf*tmask
      
!!!$omp parallel default(none) private(mytid,jj,ji) shared(jk,jpk,jpj,jpi,sdta,tmask,buf)
#ifdef __OPENMP1
        mytid = omp_get_thread_num()  ! take the thread ID
#endif
         IF (jk+mytid <=jpk) then
            DO jj= 1,jpj
            DO ji= 1,jpi
              sdta(jk,jj,ji+mytid,2) = buf(jk,jj,ji+mytid)*tmask(jk,jj,ji+mytid)
            ENDDO
            ENDDO
         ENDIF
!!!$omp end parallel
      END DO

      call readnc_slice_float(nomefile,'e3t',buf)
      
      DO jk=1,jpk,ntids!sdta(:,:,:,2) = buf*tmask
      
!!!$omp parallel default(none) private(mytid,jj,ji) shared(jk,jpk,jpj,jpi,e3tdta,tmask,buf)
#ifdef __OPENMP1
        mytid = omp_get_thread_num()  ! take the thread ID
#endif
         IF (jk+mytid <=jpk) then
            DO jj= 1,jpj
            DO ji= 1,jpi
              e3tdta(jk,jj,ji+mytid,2) = buf(jk,jj,ji+mytid)*tmask(jk,jj,ji+mytid)
            ENDDO
            ENDDO
         ENDIF
!!!$omp end parallel
      END DO



      call readnc_slice_float_2d(nomefile,'sowindsp',buf2)
       flxdta(:,:,jpwind,2) = buf2*tmask(:,:,1)
      
      call readnc_slice_float_2d(nomefile,'soshfldo',buf2)
       flxdta(:,:,jpqsr ,2) = buf2*tmask(:,:,1)
      
                                                            flxdta(:,:,jpice ,2) = 0.
      call EXISTVAR(nomefile,'sowaflcd',B)
      if (B) then
      call readnc_slice_float_2d(nomefile,'sowaflcd',buf2)
       flxdta(:,:,jpemp ,2) = buf2*tmask(:,:,1)
      
      else
         if(lwp) write(*,*) 'Evaporation data not found. Forced to zero.'
         flxdta(:,:,jpemp ,2) = 0.
      endif
!!!!!!!!!!!!!!!!!

      call EXISTVAR(nomefile,'sossheiu',B)
      if (B) then
         call readnc_slice_float_2d(nomefile,'sossheiu',buf2)
         DO jj=1,jpj
           DO ji=1,jpi
            if (umask(jj,ji,1) .EQ. 1.)  flxdta(jj,ji,8 ,2) = buf2(jj,ji)
      
           END DO
         END DO
      e3u(:,:,1) = flxdta(:,:,8 ,2)
      
      else
!     Do nothing leave the init value --> domrea
      endif
!!!!!!!!!!!!!!!!!
! epascolo warning
      call EXISTVAR(nomefile,'sossheiv',B)
      if (B) then
         call readnc_slice_float_2d(nomefile,'sossheiv',buf2)
         DO jj=1,jpj
            DO ji=1,jpi
               if (vmask(jj,ji,1) .EQ. 1.)  flxdta(jj,ji,9 ,2) = buf2(jj,ji)
      
            END DO
         END DO
      e3v(:,:,1) = flxdta(:,:,9 ,2)
      
      else
!     Do nothing leave the init value --> domrea
      endif
!!!!!!!!!!!!!!!!!

      call EXISTVAR(nomefile,'sossheit',B)
      if (B) then
         call readnc_slice_float_2d(nomefile,'sossheit',buf2)
         DO jj=1,jpj
            DO ji=1,jpi
               if (tmask(jj,ji,1) .EQ. 1.)  flxdta(jj,ji,10 ,2) = buf2(jj,ji)
      
            END DO
         END DO
      e3t(:,:,1) = flxdta(:,:,10 ,2)
      
      else
!     Do nothing leave the init value --> domrea
      endif
!!!!!!!!!!!!!!!!!


!     CALL div()               ! Horizontal divergence
!     CALL wzv()               ! vertical velocity

!        could be written for OpenMP
          DO jk=1,jpk
            DO jj=1,jpj
              DO ji=1,jpi
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
         REAL(8) zweigh, Umzweigh

         INTEGER jk,jj,ji
         INTEGER uk, uj      ! aux variables for OpenMP

      INTEGER :: mytid, ntids! omp variables

#ifdef __OPENMP1
      INTEGER ::  omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
      EXTERNAL :: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
#endif


#ifdef __OPENMP1
      ntids = omp_get_max_threads() ! take the number of threads
      mytid = -1000000

#else
      ntids = 1
      mytid = 0
#endif

      Umzweigh  = 1.0 - zweigh

          DO jk=1,jpk, ntids
!!!$omp parallel default(none) private(mytid,jj,ji,uk)
!!!$omp&                       shared(jpk,jpj,jpi,jk,ub,un,udta, vb,vn,vdta,wn,wdta,avt,avtdta,tn,tdta,sn,sdta,
!!!$omp&                       zweigh,Umzweigh,tmask,umask,vmask,e3u,e3udta,e3v,e3vdta,e3t,e3tdta,e3w,e3wdta,e3t_back)
#ifdef __OPENMP1
         mytid = omp_get_thread_num()  ! take the thread ID
#endif
        if (mytid+jk.le.jpk) then

            uk = mytid+jk
            DO jj=1,jpj
              DO ji=1,jpi

                if (umask(uk,jj,ji) .NE. 0.0) then
                un(uk,jj,ji)  = (Umzweigh*  udta(uk,jj,ji,1) + zweigh*  udta(uk,jj,ji,2))
                e3u(uk,jj,ji) = (Umzweigh*  e3udta(uk,jj,ji,1) + zweigh*  e3udta(uk,jj,ji,2))
                endif

                if (vmask(uk,jj,ji) .NE. 0.0) then
                vn(uk,jj,ji)  = (Umzweigh*  vdta(uk,jj,ji,1) + zweigh*  vdta(uk,jj,ji,2))
                e3v(uk,jj,ji) = (Umzweigh*  e3vdta(uk,jj,ji,1) + zweigh*  e3vdta(uk,jj,ji,2))
                endif

                if (tmask(uk,jj,ji) .NE. 0.0) then

                 wn(uk,jj,ji) = (Umzweigh*  wdta(uk,jj,ji,1) + zweigh*  wdta(uk,jj,ji,2))
                avt(uk,jj,ji) = (Umzweigh*avtdta(uk,jj,ji,1) + zweigh*avtdta(uk,jj,ji,2))
                e3w(uk,jj,ji) = (Umzweigh*  e3wdta(uk,jj,ji,1) + zweigh*  e3wdta(uk,jj,ji,2))
       
                 tn(uk,jj,ji) = (Umzweigh*  tdta(uk,jj,ji,1) + zweigh*  tdta(uk,jj,ji,2))
                 sn(uk,jj,ji) = (Umzweigh*  sdta(uk,jj,ji,1) + zweigh*  sdta(uk,jj,ji,2))
                e3t_back(uk,jj,ji) = e3t(uk,jj,ji)
                e3t(uk,jj,ji) = (Umzweigh*  e3tdta(uk,jj,ji,1) + zweigh*  e3tdta(uk,jj,ji,2))
                endif ! tmask


              END DO
            END DO
      endif
!!!$omp end parallel
          END DO

          DO jj=1,jpj,ntids
!!!$omp parallel default(none) private(mytid,jk,uj,ji)
!!!$omp&                       shared(jpk,jpj,jpi,jj,flx,flxdta,
!!!$omp&                              vatm,freeze,emp,qsr,jpwind,jpice,jpemp,jpqsr,zweigh, Umzweigh,jpflx)
#ifdef __OPENMP1
         mytid = omp_get_thread_num()  ! take the thread ID
#endif

        if (mytid+jj.le.jpj) then
           uj = jj+mytid

            DO jk=1,jpflx
              DO ji=1,jpi
                flx(jk,uj,ji) = ( Umzweigh * flxdta(jk,uj,ji,1)+ zweigh     * flxdta(jk,uj,ji,2) )
              END DO
            END DO

            DO ji=1,jpi
                  vatm(uj,ji)   = flx(uj,ji,jpwind)
                  freeze(uj,ji) = flx(uj,ji,jpice)
                  emp(uj,ji)    = flx(uj,ji,jpemp)
                  qsr(uj,ji)    = flx(uj,ji,jpqsr)
!                 e3u(uj,ji,1)  = flx(uj,ji,8)
!                 e3v(uj,ji,1)  = flx(uj,ji,9)
!                 e3t(uj,ji,1)  = flx(uj,ji,10)
            END DO

        endif
!!!$omp end parallel
       END DO


      END SUBROUTINE ACTUALIZE_PHYS



! *************************************************************
!      SUBROUTINE SWAP
! *    copies index 2 in index 1
! *************************************************************

      SUBROUTINE swap_PHYS
         USE myalloc
         IMPLICIT NONE
         INTEGER jk,jj,ji,jdepth
         INTEGER :: mytid, ntids! omp variables

#ifdef __OPENMP1
      INTEGER ::  omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
      EXTERNAL :: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
#endif


#ifdef __OPENMP1
      ntids = omp_get_max_threads() ! take the number of threads
      mytid = -1000000

#else
      ntids = 1
      mytid = 0
#endif

          DO jk=1,jpk,ntids
!!!$omp parallel default(None) private(mytid,jdepth,jj,ji) shared(jpk,jpj,jpi,jk,udta,vdta,wdta,avtdta,tdta,sdta)
!!!$omp&           shared(e3u,e3udta,e3v,e3vdta,e3t,e3tdta,e3w,e3wdta)
#ifdef __OPENMP1
         mytid = omp_get_thread_num()  ! take the thread ID
#endif
         jdepth=jk+mytid
         if (jdepth <= jpk) then

            DO jj=1,jpj
              DO ji=1,jpi
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
          ENDIF
!!!$omp end parallel

          END DO

          DO jk=1,jpflx,ntids
!!!$omp parallel default(None) private(mytid,jj,ji) shared(jk,jpi,jpj,jpflx,flxdta)
#ifdef __OPENMP1
         mytid = omp_get_thread_num()  ! take the thread ID
#endif
           if (jk+mytid.le.jpflx) then
            DO ji=1,jpi
              DO jj=1,jpj
                flxdta(jk,jj,ji+mytid,1) = flxdta(jk,jj,ji+mytid,2)
              END DO
            END DO
           endif
!!!$omp end parallel
          END DO



      END SUBROUTINE swap_PHYS
