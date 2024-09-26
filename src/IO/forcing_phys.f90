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
      integer  :: jk,jj,ji, jstart
      ! LOCAL
      character(LEN=30) nomefile
      character(LEN=36) DeltaT_name
      double precision ssh(jpj,jpi)
      double precision diff_e3t(jpk,jpj,jpi)
      double precision, dimension(jpj,jpi)   :: e1u_x_e2u, e1v_x_e2v, e1t_x_e2t
      double precision correction_e3t, s0,s1,s2

      if (variable_rdt) then
          DeltaT_name="DELTA_T/DeltaT_"//datestring//".txt"
          open(3333,file=DeltaT_name, form="formatted")
          read(3333,'(I5)') imposed_deltaT(2)
          close(3333)
          if (lwp) write(*,*) 'Delta T = ', imposed_deltaT(2), 'seconds'
          jk = minval(imposed_deltaT)
          rdt = real(jk , 8)
      endif
      nomefile='FORCINGS/U19951206-12:00:00.nc'

! Starting I/O
! U  *********************************************************
      nomefile = 'FORCINGS/U'//datestring//'.nc'
      if(lwp) write(*,'(A,I4,A,A)') "LOAD_PHYS --> I am ", myrank, " starting reading forcing fields from ", nomefile(1:30)
      call readnc_slice_float(nomefile,'vozocrtx',buf,ingv_lon_shift)
      udta(:,:,:,2) = buf * umask
      !$acc update device(udta(:,:,:,2)) async(1)

! V *********************************************************
      nomefile = 'FORCINGS/V'//datestring//'.nc'
      call readnc_slice_float(nomefile,'vomecrty',buf,ingv_lon_shift)
      vdta(:,:,:,2) = buf * vmask
      !$acc update device(vdta(:,:,:,2)) async(1)

! W *********************************************************


      nomefile = 'FORCINGS/W'//datestring//'.nc'

      call readnc_slice_float(nomefile,'votkeavt',buf,ingv_lon_shift)
      avtdta(:,:,:,2) = buf*tmask
      !$acc update device(avtdta(:,:,:,2)) async(1)

! T *********************************************************
      nomefile = 'FORCINGS/T'//datestring//'.nc'
      call readnc_slice_float(nomefile,'votemper',buf,ingv_lon_shift)
      tdta(:,:,:,2) = buf*tmask
      !$acc update device(tdta(:,:,:,2)) async(1)

      call readnc_slice_float(nomefile,'vosaline',buf,ingv_lon_shift)
      sdta(:,:,:,2) = buf*tmask
      !$acc update device(sdta(:,:,:,2)) async(1)

    if (IS_FREE_SURFACE) then
         call readnc_slice_float_2d(nomefile,'sossheig',buf2,ingv_lon_shift)
         ssh = buf2*tmask(1,:,:)

          e3tdta(:,:,:,2) = e3t_0
          DO ji= 1,jpi
          DO jj= 1,jpj
          if (tmask(1,jj,ji).eq.1) then  ! to do the division
              correction_e3t=( 1.0 + ssh(jj,ji)/h_column(jj,ji))
              DO jk=1,mbathy(jj,ji)
                   e3tdta(jk,jj,ji,2)  = e3t_0(jk,jj,ji) * correction_e3t
              ENDDO
          endif
          ENDDO
          ENDDO

         e1u_x_e2u = e1u*e2u
         e1v_x_e2v = e1v*e2v
         e1t_x_e2t = e1t*e2t

         diff_e3t = e3tdta(:,:,:,2) - e3t_0
         e3udta(:,:,:,2) = 0.0
         e3vdta(:,:,:,2) = 0.0

         DO ji = 1,jpim1
         DO jj = 1,jpjm1
         DO jk = 1,jpk
             s0= e1t_x_e2t(jj,ji ) * diff_e3t(jk,jj,ji)
             s1= e1t_x_e2t(jj,ji+1) * diff_e3t(jk,jj,ji+1)
             s2= e1t_x_e2t(jj+1,ji) * diff_e3t(jk,jj+1,ji)
             e3udta(jk,jj,ji,2) = 0.5*(umask(jk,jj,ji)/(e1u_x_e2u(jj,ji)) * (s0 + s1))
             e3vdta(jk,jj,ji,2) = 0.5*(vmask(jk,jj,ji)/(e1v_x_e2v(jj,ji)) * (s0 + s2))
         ENDDO
         ENDDO
         ENDDO

         DO ji = 1,jpi
         DO jj = 1,jpj
         DO jk = 1,jpk
             e3udta(jk,jj,ji,2) = e3u_0(jk,jj,ji) + e3udta(jk,jj,ji,2)
             e3vdta(jk,jj,ji,2) = e3v_0(jk,jj,ji) + e3vdta(jk,jj,ji,2)
         ENDDO
         ENDDO
         ENDDO



         DO ji = 1,jpi
         DO jj = 1,jpj
             e3wdta(1,jj,ji,2) = e3w_0(1,jj,ji) + diff_e3t(1,jj,ji)
         ENDDO
         ENDDO

         DO ji = 1,jpi
         DO jj = 1,jpj
         DO jk = 2,mbathy(jj,ji)
              e3wdta(jk,jj,ji,2) = e3w_0(jk,jj,ji) + 0.5*( diff_e3t(jk-1,jj,ji) + diff_e3t(jk,jj,ji))
         ENDDO
         jstart = jk
         DO jk =  jstart, jpk
             e3wdta(jk,jj,ji,2) = e3w_0(jk,jj,ji) + diff_e3t(jk-1,jj,ji)
         ENDDO

         ENDDO
         ENDDO

         !$acc update device(e3udta(:,:,:,2),e3vdta(:,:,:,2),e3wdta(:,:,:,2)) async(1)
         !$acc update device(e3tdta(:,:,:,2)) async(1)

     endif ! IS_FREE_SURFACE




      if (ingv_files_direct_reading) then
           nomefile = 'FORCINGS/U'//datestring//'.nc'
           call readnc_slice_float_2d(nomefile,'sozotaux',buf2,ingv_lon_shift)
           taux = buf2*tmask(1,:,:)*umask(1,:,:)

           nomefile = 'FORCINGS/V'//datestring//'.nc'
           call readnc_slice_float_2d(nomefile,'sometauy',buf2,ingv_lon_shift)
           tauy = buf2*tmask(1,:,:)*vmask(1,:,:)

           call PURE_WIND_SPEED(taux,tauy,jpi,jpj, buf2)
      else
          nomefile = 'FORCINGS/T'//datestring//'.nc'
          call readnc_slice_float_2d(nomefile,'sowindsp',buf2,ingv_lon_shift)
      endif
      flxdta(:,:,jpwind,2) = buf2*tmask(1,:,:) * spongeT

      nomefile = 'FORCINGS/T'//datestring//'.nc'
      call readnc_slice_float_2d(nomefile,'soshfldo',buf2,ingv_lon_shift)
      flxdta(:,:,jpqsr ,2) = buf2*tmask(1,:,:) * spongeT
      flxdta(:,:,jpice ,2) = 0.
      flxdta(:,:,jpemp ,2) = 0.


      if (read_W_from_file) then
          nomefile = 'FORCINGS/W'//datestring//'.nc'
          call readnc_slice_float(nomefile,'vovecrtz',buf,ingv_lon_shift)
          wdta(:,:,:,2) = buf * tmask
      else
          CALL COMPUTE_W()               ! vertical velocity
      endif
      !$acc update device(wdta(:,:,:,2)) async(1)
      
!        could be written for OpenMP
      !$acc parallel loop gang vector default(present) collapse(3) async(1)
      DO ji=1,jpi
         DO jj=1,jpj
            DO jk=1,jpk
               udta(jk,jj,ji,2) =   udta(jk,jj,ji,2) * spongeVel(jk,jj,ji)
               vdta(jk,jj,ji,2) =   vdta(jk,jj,ji,2) * spongeVel(jk,jj,ji)
               wdta(jk,jj,ji,2) =   wdta(jk,jj,ji,2) * spongeVel(jk,jj,ji)
               avtdta(jk,jj,ji,2) = avtdta(jk,jj,ji,2) * spongeVel(jk,jj,ji)
               tn(jk,jj,ji)=tdta(jk,jj,ji,2)
               sn(jk,jj,ji)=sdta(jk,jj,ji,2)
            END DO
         END DO
      END DO
      !$acc end parallel
      !$acc wait(1)

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


         !$acc parallel loop gang vector default(present) collapse(3) async(1)
         DO ji=1,jpi
            DO jj=1,jpj
               DO jk=1,jpk
                  un(jk,jj,ji) = Umzweigh* udta(jk,jj,ji,1) + zweigh*  udta(jk,jj,ji,2)
                  vn(jk,jj,ji) = Umzweigh* vdta(jk,jj,ji,1) + zweigh*  vdta(jk,jj,ji,2)
                  wn(jk,jj,ji) = Umzweigh* wdta(jk,jj,ji,1) + zweigh*  wdta(jk,jj,ji,2)
                  tn(jk,jj,ji)  = Umzweigh* tdta(jk,jj,ji,1)   + zweigh* tdta(jk,jj,ji,2)
                  sn(jk,jj,ji)  = Umzweigh* sdta(jk,jj,ji,1)   + zweigh* sdta(jk,jj,ji,2)
                  avt(jk,jj,ji) = Umzweigh* avtdta(jk,jj,ji,1) + zweigh* avtdta(jk,jj,ji,2)
               ENDDO
            ENDDO
         ENDDO
         !$acc end parallel

         IF (IS_FREE_SURFACE) then
            !$acc parallel loop gang vector default(present) collapse(3) async(1)
            DO ji=1,jpi
               DO jj=1,jpj
                  DO jk=1,jpk
                     e3u(jk,jj,ji) = Umzweigh* e3udta(jk,jj,ji,1) + zweigh* e3udta(jk,jj,ji,2)
                     e3v(jk,jj,ji) = Umzweigh* e3vdta(jk,jj,ji,1) + zweigh* e3vdta(jk,jj,ji,2)
                     e3w(jk,jj,ji) = Umzweigh* e3wdta(jk,jj,ji,1) + zweigh* e3wdta(jk,jj,ji,2)
                  ENDDO
               ENDDO
            ENDDO
            !$acc end parallel

            if (forcing_phys_initialized) then
               !$acc update device(e3t) if(forcing_phys_initialized) async(1)
               !$acc parallel loop gang vector default(present) collapse(3) async(1)
               DO ji=1,jpi
                  DO jj=1,jpj
                     DO jk=1,jpk
                        e3t_back(jk,jj,ji) = e3t(jk,jj,ji)
                        e3t(jk,jj,ji) = (Umzweigh*  e3tdta(jk,jj,ji,1) + zweigh*  e3tdta(jk,jj,ji,2))
                     ENDDO
                  ENDDO
               ENDDO
               !$acc end parallel
            else
               !$acc parallel loop gang vector default(present) collapse(3) async(1)
               DO ji=1,jpi
                  DO jj=1,jpj
                     DO jk=1,jpk
                        e3t(jk,jj,ji) = (Umzweigh*  e3tdta(jk,jj,ji,1) + zweigh*  e3tdta(jk,jj,ji,2))
                        e3t_back(jk,jj,ji) = e3t(jk,jj,ji)
                     ENDDO
                  ENDDO
               ENDDO
               !$acc end parallel
               forcing_phys_initialized = .TRUE.
            endif
        ENDIF

        flx = Umzweigh * flxdta(:,:,:,1) + zweigh * flxdta(:,:,:,2)

        !$acc update device(flx) async(1)
        !$acc parallel loop gang vector default(present) collapse(2) async(1)
        DO ji=1,jpi
           DO uj=1,jpj
              vatm(uj,ji)   = flx(uj,ji,jpwind)
              emp(uj,ji)    = flx(uj,ji,jpemp)
              qsr(uj,ji)    = flx(uj,ji,jpqsr)
           END DO
        END DO
        !$acc end parallel loop

        DO ji=1,jpi
           DO uj=1,jpj
              freeze(uj,ji) = flx(uj,ji,jpice)
           END DO
        END DO

        !$acc wait(1)

      END SUBROUTINE ACTUALIZE_PHYS



! *************************************************************
!     SUBROUTINE SWAP
!     copies index 2 in index 1
! **************************************************************

      SUBROUTINE swap_PHYS
         USE myalloc
         IMPLICIT NONE
         integer :: ji,jj,jk

         !$acc parallel loop gang vector default(present) collapse(3)
         DO ji=1,jpi
            DO jj=1,jpj
               DO jk=1,jpk
                  udta(jk,jj,ji,1) =    udta(jk,jj,ji,2)
                  vdta(jk,jj,ji,1) =    vdta(jk,jj,ji,2)
                  wdta(jk,jj,ji,1) =    wdta(jk,jj,ji,2)
                  avtdta(jk,jj,ji,1) =  avtdta(jk,jj,ji,2)
                  tdta(jk,jj,ji,1) =    tdta(jk,jj,ji,2)
                  sdta(jk,jj,ji,1) =    sdta(jk,jj,ji,2)
               ENDDO
            ENDDO
         ENDDO
         !$acc end parallel
         flxdta(:,:,:,1) =  flxdta(:,:,:,2)

         IF (IS_FREE_SURFACE) then
            !$acc parallel loop gang vector default(present) collapse(3)
            DO ji=1,jpi
               DO jj=1,jpj
                  DO jk=1,jpk
                     e3udta(jk,jj,ji,1) =  e3udta(jk,jj,ji,2)
                     e3vdta(jk,jj,ji,1) =  e3vdta(jk,jj,ji,2)
                     e3wdta(jk,jj,ji,1) =  e3wdta(jk,jj,ji,2)
                     e3tdta(jk,jj,ji,1) =  e3tdta(jk,jj,ji,2)
                  ENDDO
               ENDDO
            ENDDO
            !$acc end parallel
         ENDIF

         imposed_deltaT(1) = imposed_deltaT(2)
      END SUBROUTINE swap_PHYS

! ************************************************
!     INIT_PHYS
!     prepares sponge variables
! ************************************************
      SUBROUTINE INIT_PHYS
      USE myalloc

! ----------------------------------------------------------------------
!  BEGIN BC_REFACTORING SECTION
!  ---------------------------------------------------------------------

      use BC_mem
      use bc_set_mod

! ----------------------------------------------------------------------
!  END BC_REFACTORING SECTION
!  ---------------------------------------------------------------------

      IMPLICIT NONE
      ! integer ji, jj,jk
      ! double precision reduction_value, alpha
      ! double precision lon_limit

      ! lon_limit = -7.5
      ! alpha     = 1.0
      spongeT     = 1.0
      spongeVel   = 1.0

      ! if (internal_sponging) then
      !     DO ji=1,jpi
      !     DO jj=1,jpj
      !         if (glamt(jj,ji).lt.lon_limit) then
      !             reduction_value = 1.e-6
      !             spongeT(jj,ji) = reduction_value
      !         endif
      !     ENDDO
      !     ENDDO


      !     DO ji=1,jpi
      !     DO jj=1,jpj
      !         if (glamt(jj,ji).lt.lon_limit) then
      !             reduction_value = exp( -alpha*(  (glamt(jj,ji)-lon_limit)**2)  )
      !             do jk=1,jpk
      !                 spongeVel(jk,jj,ji) = reduction_value
      !             enddo
      !         endif
      !     ENDDO
      !     ENDDO
      ! endif

! ----------------------------------------------------------------------
!  BEGIN BC_REFACTORING SECTION
!  ---------------------------------------------------------------------

      call boundaries%apply_phys(glamt, spongeT, spongeVel, internal_sponging)
      !$acc update device(spongeVel)

! ----------------------------------------------------------------------
!  END BC_REFACTORING SECTION
!  ---------------------------------------------------------------------


      END SUBROUTINE INIT_PHYS


      SUBROUTINE COMPUTE_W ()
!---------------------------------------------------------------------
!
!                       ROUTINE compute_w
!                     ***************
!
!  Purpose :
!  ---------
!   Compute the now vertical velocity after the array swap.
!
!   Method :
!   -------
!   Using the incompressibility hypothesis, the vertical velocity
!   is computed by integrating the horizontal divergence from the
!   bottom to the surface.
!   The boundary conditions are w=0 at the bottom (no flux) and
!   w=0 at the sea surface (rigid lid).
!

       USE myalloc
       USE ogstm_mpi_module
       IMPLICIT NONE

!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER ji, jj, jk


      double precision zbt
      double precision zwu(jpj,jpi), zwv(jpj,jpi)
      double precision hdivn(jpk,jpj,jpi)

      hdivn = 0.0


      IF (IS_FREE_SURFACE) THEN

      DO jk = 1, jpkm1


! 1. Horizontal fluxes

        DO jj = 1, jpjm1
        DO ji = 1, jpim1
            zwu(jj,ji) = e2u(jj,ji) * e3udta(jk,jj,ji,2) * udta(jk,jj,ji,2)
            zwv(jj,ji) = e1v(jj,ji) * e3vdta(jk,jj,ji,2) * vdta(jk,jj,ji,2)
        END DO
        END DO


! 2. horizontal divergence


        DO jj = 2, jpjm1
        DO ji = 2, jpim1
            zbt = e1t(jj,ji) * e2t(jj,ji) * e3tdta(jk,jj,ji,2)
            hdivn(jk,jj,ji) = (  zwu(jj,ji) - zwu(jj,ji-1  ) &
                               + zwv(jj,ji) - zwv(jj-1  ,ji)  ) / zbt
        END DO
        END DO


      END DO
      ELSE
      DO jk = 1, jpkm1


! 1. Horizontal fluxes

        DO jj = 1, jpjm1
        DO ji = 1, jpim1
            zwu(jj,ji) = e2u(jj,ji) * e3u(jk,jj,ji) * udta(jk,jj,ji,2)
            zwv(jj,ji) = e1v(jj,ji) * e3v(jk,jj,ji) * vdta(jk,jj,ji,2)
        END DO
        END DO


! 2. horizontal divergence


        DO jj = 2, jpjm1
        DO ji = 2, jpim1
            zbt = e1t(jj,ji) * e2t(jj,ji) * e3t(jk,jj,ji)
            hdivn(jk,jj,ji) = (  zwu(jj,ji) - zwu(jj,ji-1  ) &
                               + zwv(jj,ji) - zwv(jj-1  ,ji)  ) / zbt
        END DO
        END DO


      END DO

      ENDIF


! 3. Lateral boundary conditions on hdivn

#ifdef key_mpp
! ... Mpp : export boundary values to neighboring processors
     CALL mpplnk_my( hdivn )
#endif




! 1. Surface and bottom boundary condition: w=0 (rigid lid and no flux)
! ----------------------------------------
      wdta(  1,:,:,2 ) = 0.0
      wdta(jpk,:,:,2 ) = 0.0

! 2. Computation from the bottom
! ------------------------------
IF (IS_FREE_SURFACE) THEN
     DO ji = 1, jpi
     DO jj = 1, jpj
     DO jk = jpkm1, 1, -1

            wdta(jk,jj,ji,2) = wdta(jk+1,jj,ji,2) - e3tdta(jk,jj,ji,2)*hdivn(jk,jj,ji)

     END DO
     END DO
     END DO
ELSE
     DO ji = 1, jpi
     DO jj = 1, jpj
     DO jk = jpkm1, 1, -1

            wdta(jk,jj,ji,2) = wdta(jk+1,jj,ji,2) - e3t(jk,jj,ji)*hdivn(jk,jj,ji)

     END DO
     END DO
     END DO
ENDIF

END SUBROUTINE COMPUTE_W


SUBROUTINE PURE_WIND_SPEED(TAUx, TAUy, sizex,sizey, WSP)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: sizex, sizey
    double precision, DIMENSION(sizex,sizey), INTENT(IN ) :: TAUx, TAUy
    double precision, DIMENSION(sizex,sizey), INTENT(OUT) :: WSP
    ! local
    double precision rho, Cdrag, K
    rho    = 1.3 ! kg/m3
    Cdrag  = 1.5 * 0.001

    K     = sqrt(1/(rho*Cdrag));
    WSP = (TAUx**2 + TAUy**2)**0.25 * K


END SUBROUTINE PURE_WIND_SPEED

