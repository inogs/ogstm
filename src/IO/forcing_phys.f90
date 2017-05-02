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
      LOGICAL :: B, IS_INGV_E3T
      integer  :: jk,jj,ji, jstart
      ! LOCAL
      character(LEN=30) nomefile
      double precision ssh(jpj,jpi)
      double precision diff_e3t(jpk,jpj,jpi)
      double precision, dimension(jpj,jpi)   :: e1u_x_e2u, e1v_x_e2v, e1t_x_e2t
      double precision correction_e3t, s0,s1,s2

      nomefile='FORCINGS/U19951206-12:00:00.nc'

! Starting I/O
! U  *********************************************************
      nomefile = 'FORCINGS/U'//datestring//'.nc'
      if(lwp) write(*,'(A,I4,A,A)') "LOAD_PHYS --> I am ", myrank, " starting reading forcing fields from ", nomefile(1:30)
      call readnc_slice_float(nomefile,'vozocrtx',buf,ingv_lon_shift)
      udta(:,:,:,2) = buf * umask * spongeVel

      IS_INGV_E3T = .false. ! call EXISTVAR(nomefile,'e3u',IS_INGV_E3T)
      if (IS_INGV_E3T) then
          call readnc_slice_float(nomefile,'e3u',buf,ingv_lon_shift)
          e3udta(:,:,:,2) = buf!*umask
      endif



! V *********************************************************
      nomefile = 'FORCINGS/V'//datestring//'.nc'
      call readnc_slice_float(nomefile,'vomecrty',buf,ingv_lon_shift)
      vdta(:,:,:,2) = buf*vmask * spongeVel
      

      if (IS_INGV_E3T) then
          call readnc_slice_float(nomefile,'e3v',buf,ingv_lon_shift)
          e3vdta(:,:,:,2) = buf!*vmask
      endif



! W *********************************************************


      nomefile = 'FORCINGS/W'//datestring//'.nc'

!      call readnc_slice_float(nomefile,'vovecrtz',buf)
!      wdta(:,:,:,2) = buf * tmask * spongeVel

      call readnc_slice_float(nomefile,'votkeavt',buf,ingv_lon_shift)
      avtdta(:,:,:,2) = buf*tmask

      if (IS_INGV_E3T) then
          call readnc_slice_float(nomefile,'e3w',buf,ingv_lon_shift)
          e3wdta(:,:,:,2) = buf!*tmask
      endif


! T *********************************************************
      nomefile = 'FORCINGS/T'//datestring//'.nc'
      call readnc_slice_float(nomefile,'votemper',buf,ingv_lon_shift)
      tdta(:,:,:,2) = buf*tmask

      call readnc_slice_float(nomefile,'vosaline',buf,ingv_lon_shift)
      sdta(:,:,:,2) = buf*tmask


      if (IS_INGV_E3T) then
          call readnc_slice_float(nomefile,'e3t',buf,ingv_lon_shift)
          e3tdta(:,:,:,2) = buf!*tmask
      endif

    if (.not.IS_INGV_E3T) then
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




     endif ! IS_INGV_E3T




      if (ingv_files_direct_reading) then
           nomefile = 'FORCINGS/U'//datestring//'.nc'
           call readnc_slice_float_2d(nomefile,'sozotaux',buf2,ingv_lon_shift)
           taux = buf2*tmask(1,:,:)

           nomefile = 'FORCINGS/V'//datestring//'.nc'
           call readnc_slice_float_2d(nomefile,'sometauy',buf2,ingv_lon_shift)
           tauy = buf2*tmask(1,:,:)

           call PURE_WIND_SPEED(taux,tauy,jpi,jpj, buf2)
      else
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
          wdta(:,:,:,2) = buf * tmask * spongeVel
      else
          CALL COMPUTE_W()               ! vertical velocity
      endif



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

         un = Umzweigh* udta(:,:,:,1) + zweigh*  udta(:,:,:,2)
         vn = Umzweigh* vdta(:,:,:,1) + zweigh*  vdta(:,:,:,2)
         wn = Umzweigh* wdta(:,:,:,1) + zweigh*  wdta(:,:,:,2)

         e3u = Umzweigh* e3udta(:,:,:,1) + zweigh* e3udta(:,:,:,2)
         e3v = Umzweigh* e3vdta(:,:,:,1) + zweigh* e3vdta(:,:,:,2)
         e3w = Umzweigh* e3wdta(:,:,:,1) + zweigh* e3wdta(:,:,:,2)

         tn  = Umzweigh* tdta(:,:,:,1)   + zweigh* tdta(:,:,:,2)
         sn  = Umzweigh* sdta(:,:,:,1)   + zweigh* sdta(:,:,:,2)
         avt = Umzweigh* avtdta(:,:,:,1) + zweigh* avtdta(:,:,:,2)

        if (forcing_phys_initialized) then
           e3t_back = e3t
           e3t = (Umzweigh*  e3tdta(:,:,:,1) + zweigh*  e3tdta(:,:,:,2))
        else
          e3t = (Umzweigh*  e3tdta(:,:,:,1) + zweigh*  e3tdta(:,:,:,2))
          e3t_back = e3t
          forcing_phys_initialized = .TRUE.
        endif

        flx = Umzweigh * flxdta(:,:,:,1) + zweigh * flxdta(:,:,:,2)



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

                    udta(:,:,:,1) =    udta(:,:,:,2)
                  e3udta(:,:,:,1) =  e3udta(:,:,:,2)
                    vdta(:,:,:,1) =    vdta(:,:,:,2)
                  e3vdta(:,:,:,1) =  e3vdta(:,:,:,2)
                    wdta(:,:,:,1) =    wdta(:,:,:,2)
                  e3wdta(:,:,:,1) =  e3wdta(:,:,:,2)
                  avtdta(:,:,:,1) =  avtdta(:,:,:,2)
                    tdta(:,:,:,1) =    tdta(:,:,:,2)
                    sdta(:,:,:,1) =    sdta(:,:,:,2)
                  e3tdta(:,:,:,1) =  e3tdta(:,:,:,2)
                  flxdta(:,:,:,1) =  flxdta(:,:,:,2)


      END SUBROUTINE swap_PHYS

! ************************************************
!     INIT_PHYS
!     prepares sponge variables
! ************************************************
      SUBROUTINE INIT_PHYS
      USE myalloc
      IMPLICIT NONE
      integer ji, jj,jk
      double precision reduction_value, alpha
      double precision lon_limit

      lon_limit = -7.5
      alpha     = 1.0
      spongeT     = 1.0
      spongeVel   = 1.0

      if (internal_sponging) then
          DO ji=1,jpi
          DO jj=1,jpj
              if (glamt(jj,ji).lt.lon_limit) then
                  reduction_value = 1.e-6
                  spongeT(jj,ji) = reduction_value
              endif
          ENDDO
          ENDDO


          DO ji=1,jpi
          DO jj=1,jpj
              if (glamt(jj,ji).lt.lon_limit) then
                  reduction_value = exp( -alpha*(  (glamt(jj,ji)-lon_limit)**2)  )
                  do jk=1,jpk
                      spongeVel(jk,jj,ji) = reduction_value
                  enddo
              endif
          ENDDO
          ENDDO
      endif


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
     DO ji = 1, jpi
     DO jj = 1, jpj
     DO jk = jpkm1, 1, -1

            wdta(jk,jj,ji,2) = wdta(jk+1,jj,ji,2) - e3tdta(jk,jj,ji,2)*hdivn(jk,jj,ji)

     END DO
     END DO
     END DO


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

