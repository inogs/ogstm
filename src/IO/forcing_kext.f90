      SUBROUTINE forcings_KEXT(datestring)
!---------------------------------------------------------------------


       USE myalloc
       ! epascolo USE myalloc_mpp
       USE OPT_mem
       USE TIME_MANAGER
       USE mpi
       IMPLICIT NONE

      character(LEN=17), INTENT(IN) ::  datestring

! local declarations
! ==================
      double precision sec,zweigh
      integer Before, After
      INTEGER iswap

      forcing_kext_partTime = MPI_WTIME()

      sec=datestring2sec(DATEstring)
      call TimeInterpolation(sec,TC_OPTATM, BEFORE, AFTER, zweigh)

      iswap  = 0


! ----------------------- INITIALIZATION -------------
      IF (datestring.eq.DATESTART) then

          CALL LOAD_KEXT(TC_OPTATM%TimeStrings(TC_OPTATM%Before)) ! CALL dynrea(iperm1)
          iswap = 1
          call swap_KEXT


        CALL LOAD_KEXT(TC_OPTATM%TimeStrings(TC_OPTATM%After)) !CALL dynrea(iper)


      ENDIF





! --------------------------------------------------------
! and now what we have to DO at every time step
! --------------------------------------------------------

! check the validity of the period in memory

      IF (BEFORE.ne.TC_OPTATM%Before) then
         TC_OPTATM%Before = BEFORE
         TC_OPTATM%After  = AFTER

         call swap_KEXT
         iswap = 1


          CALL LOAD_KEXT(TC_OPTATM%TimeStrings(TC_OPTATM%After))

          IF(lwp) WRITE (numout,*) ' Extinction factor DATA READ for Time = ', TC_OPTATM%TimeStrings(TC_OPTATM%After)
!      ******* LOADED NEW FRAME *************
      END IF





! compute the DATA at the given time step

      SELECT CASE (nsptint)
           CASE (0)  !  ------- no time interpolation
!      we have to initialize DATA IF we have changed the period
              IF (iswap.eq.1) THEN
                 zweigh = 1.0
                 call ACTUALIZE_KEXT(zweigh)! initialize now fields with the NEW DATA READ
              END IF

          CASE (1) ! ------------linear interpolation ---------------
             call ACTUALIZE_KEXT(zweigh)
      END SELECT


      forcing_kext_partTime = MPI_WTIME() - forcing_kext_partTime
      forcing_kext_TotTime  = forcing_kext_TotTime + forcing_kext_partTime



      RETURN
      END

! ******************************************************
!     SUBROUTINE LOAD_KEXT(datestring)
!     loads OPTICS/atm_yyyy0107-00:00:00.nc
! ******************************************************
      SUBROUTINE LOAD_KEXT(datestring)
! ======================
      USE calendar
      USE myalloc
      ! epascolo USE myalloc_mpp
      USE OPT_mem
      USE TIME_MANAGER
      USE BC_mem

      IMPLICIT NONE

      CHARACTER(LEN=17), INTENT(IN) :: datestring

      character(LEN=39) nomefile,fileformat
      character(LEN=4) yyyy
      character(LEN=2) mm
      double precision ::  junk(jpj,jpi)

      nomefile='OPTICS/yyyy/mm/atm.yyyy0107-00:00:00.nc'
      fileformat='("OPTICS/",A4,"/",A2,"/atm.",A17,".nc")'


!     Starting I/O
!    **********************************************************
      yyyy=datestring(1:4)
      mm=datestring(5:6)
      WRITE ( nomefile, fileformat ) yyyy,mm,datestring

      if(lwp) write(*,'(A,I4,A,A)') "LOAD_KEXT --> I am ", myrank, " starting reading atmospheric fields from ", nomefile

       call readnc_slice_float_2d(nomefile,'sp',buf2,0)
       spIO(:,:,2) = buf2*tmask(1,:,:)

       call readnc_slice_float_2d(nomefile,'msl',buf2,0)
       mslIO(:,:,2) = buf2*tmask(1,:,:)

       call readnc_slice_float_2d(nomefile,'t2m',buf2,0)
       t2mIO(:,:,2) = buf2*tmask(1,:,:)

       call readnc_slice_float_2d(nomefile,'d2m',buf2,0)
       d2mIO(:,:,2) = buf2*tmask(1,:,:)

       call readnc_slice_float_2d(nomefile,'tcc',buf2,0)
       tccIO(:,:,2) = buf2*tmask(1,:,:)



       call readnc_slice_float_2d(nomefile,'u10',buf2,0)
       call readnc_slice_float_2d(nomefile,'v10',junk,0)
       buf2 = sqrt(buf2*buf2 + junk*junk)
       w10IO(:,:,2) = buf2*tmask(1,:,:)

       call readnc_slice_float_2d(nomefile,'tclw',buf2,0)
       tclwIO(:,:,2) = buf2*tmask(1,:,:)
       call readnc_slice_float_2d(nomefile,'tco3',buf2,0)
       tco3IO(:,:,2) = buf2*tmask(1,:,:)



      

      END SUBROUTINE LOAD_KEXT

! ******************************************************
!     SUBROUTINE ACTUALIZE_PHYS(zweigh)
!     performs time interpolation
!     x(1)*(1-zweigh) + x(2)*zweigh
! ******************************************************

    SUBROUTINE actualize(zweigh,array3d, array2d)
    use myalloc
    double precision, INTENT(IN) :: zweigh
    double precision, INTENT(IN) :: array3d(jpj,jpi,2)
    double precision, INTENT(OUT) :: array2d(jpj,jpi)

    double precision Umzweigh

    INTEGER jj,ji

    Umzweigh  = 1.0 - zweigh
        DO ji=1,jpi
          DO jj=1,jpj
             if (tmask(1,jj,ji).ne.0 ) then
             array2d(jj,ji) = ( Umzweigh * array3d(jj,ji,1)+ zweigh * array3d(jj,ji,2) )
             else
             array2d(jj,ji) = 1.0
             endif

          END DO
        END DO
   END SUBROUTINE actualize

    SUBROUTINE swap(array3d)
    use myalloc
    double precision :: array3d(jpj,jpi,2)

    INTEGER jj,ji

        DO ji=1,jpi
          DO jj=1,jpj
             array3d(jj,ji,1) = array3d(jj,ji,2)
          END DO
        END DO
   END SUBROUTINE swap


! ******************************************************
!     SUBROUTINE ACTUALIZE_PHYS(zweigh)
!     performs time interpolation
!     x(1)*(1-zweigh) + x(2)*zweigh
! ******************************************************
      SUBROUTINE ACTUALIZE_KEXT(zweigh)
         USE myalloc
         USE OPT_mem
         IMPLICIT NONE

         double precision, INTENT(IN) :: zweigh


        call actualize(zweigh,spIO,sp)
        call actualize(zweigh,mslIO,msl)
        call actualize(zweigh,t2mIO,t2m)
        call actualize(zweigh,d2mIO,d2m)
        call actualize(zweigh,tccIO,tcc)
        call actualize(zweigh,w10IO,w10)
        call actualize(zweigh,tclwIO,tclw)
        call actualize(zweigh,tco3IO,tco3)



      END SUBROUTINE ACTUALIZE_KEXT




! *************************************************************
!      SUBROUTINE SWAP_KEXT
! *    copia l'indice 2 nell'indice 1
! *************************************************************

      SUBROUTINE swap_KEXT
         USE myalloc
         USE OPT_mem
         IMPLICIT NONE

        call swap(spIO)
        call swap(mslIO)
        call swap(t2mIO)
        call swap(d2mIO)
        call swap(tccIO)
        call swap(w10IO)
        call swap(tclwIO)
        call swap(tco3IO)


      END SUBROUTINE swap_KEXT


