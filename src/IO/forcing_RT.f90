      SUBROUTINE forcings_RT(datestring)
!---------------------------------------------------------------------
!
       USE myalloc
       USE TIME_MANAGER
       use mpi
       use OPT_mem
       IMPLICIT NONE

      character(LEN=17), INTENT(IN) ::  datestring

! local declarations
! ==================
      double precision :: sec
      double precision :: day_seconds ! seconds from the day begin
      integer          :: year, month, day_RTactual



!     iswap  : indicator of swap of dynamic DATA array

       forcing_rt_partTime = MPI_WTIME()  ! cronometer-start

! ----------------------- INITIALISATION -------------
! and now what we have to DO at every time step
! --------------------------------------------------------
! ----------------------- INITIALISATION -------------
      if (datestring.eq.DATESTART) then 
          CALL LOAD_RT(datestring)
          call read_date_string(datestring, year, month, day_RTcheck, day_seconds) 
      endif

      call read_date_string(datestring, year, month, day_RTactual, day_seconds) 



! check the validity of the period in memory

      if (day_RTactual.ne.day_RTcheck) then
!      ******* LOADED NEW FRAME *************

         CALL LOAD_RT(datestring)
         day_RTcheck = day_RTactual
      END IF

      call ACTUALIZE_RT(day_seconds)! initialize now fields with the NEW DATA READ

      call CALC_QSR()


       forcing_rt_partTime = MPI_WTIME() - forcing_rt_partTime
       forcing_rt_TotTime  = forcing_rt_TotTime  + forcing_rt_partTime


      END SUBROUTINE forcings_RT

! ******************************************************
!     SUBROUTINE LOAD_PHYS(datestring)
!
!
! ******************************************************
       SUBROUTINE LOAD_RT(datestring)
! ======================
      USE calendar
      USE myalloc
      USE TIME_MANAGER
      USE OPT_mem

      IMPLICIT NONE

      character(LEN=17), INTENT(IN) :: datestring
      LOGICAL :: B
      ! LOCAL
      character(LEN=30) nomefile

      nomefile='OASIM/rad_0m20030101.nc'

! Starting I/O
! Ed_0m  *********************************************************

      nomefile = 'OASIM/rad_0m'//TRIM(datestring(1:8))//'.nc'
      if(lwp) write(*,'(A,I4,A,A)') "LOAD_RT --> I am ", myrank, " starting reading radiative transfer atm model fields from ", nomefile(1:23)
      call readnc_OASIM_float_2d(nomefile,'Ed_0m', Ed_0m_COARSE, OASIM_lon, OASIM_lat)


! Es_0m *********************************************************

      nomefile = 'OASIM/rad_0m'//TRIM(datestring(1:8))//'.nc'
      if(lwp) write(*,'(A,I4,A,A)') "LOAD_RT --> I am ", myrank, " starting reading radiative transfer atm model fields from ", nomefile(1:23)
      call readnc_OASIM_float_2d(nomefile,'Es_0m', Es_0m_COARSE, OASIM_lon, OASIM_lat)


      END SUBROUTINE LOAD_RT

! ******************************************************
      SUBROUTINE ACTUALIZE_RT(day_seconds)
         USE myalloc
         USE OPT_mem
         IMPLICIT NONE
         double precision day_seconds

         INTEGER it ! from 1 to 12
         INTEGER jj,ji
         INTEGER jj_coarse,ji_coarse

         it=int(day_seconds/(7200.)) + 1   

         do ji=1,jpi
            do jj=1,jpj
               jj_coarse=nint(gphit(jj,ji)-OASIM_lat(1,1))
               ji_coarse=nint(glamt(jj,ji)-OASIM_lon(1,1))

               Ed_0m(:,jj,ji) = Ed_0m_COARSE(:,it,jj_coarse,ji_coarse)
               Es_0m(:,jj,ji) = Es_0m_COARSE(:,it,jj_coarse,ji_coarse)
               
            enddo
         enddo


      END SUBROUTINE ACTUALIZE_RT

      SUBROUTINE CALC_QSR()
         USE myalloc
         USE OPT_mem
         IMPLICIT NONE

         INTEGER jj,ji,wl
         do ji=1,jpi
            do jj=1,jpj

              qsr(jj,ji) = 0.0

              do wl=1,nwl
                 qsr(jj,ji) = qsr(jj,ji) + Ed_0m(wl,jj,ji) + Es_0m(wl,jj,ji) 
              enddo

            enddo
         enddo

        qsr = qsr * spongeT


      END SUBROUTINE CALC_QSR


