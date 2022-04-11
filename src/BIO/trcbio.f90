
      SUBROUTINE trcbio
!!!---------------------------------------------------------------------
!!!
!!!                       ROUTINE trcbio
!!!                     *******************
!!!
!!!  PURPOSE :
!!!  ---------
!!!     compute the now trend due to biogeochemical processes
!!!     and add it to the general trend of passive tracers equations.
!!!
!!!    Three options:
!!!
!!!   METHOD :
!!!   -------
!!!      each now biological flux is calculated  in FUNCTION of now
!!!      concentrations of tracers.
!!!      depending on the tracer, these fluxes are sources or sinks.
!!!      the total of the sources and sinks for each tracer
!!!      is added to the general trend.
!!!
!!!        tra = tra + zf...tra - zftra...
!!!                             |         |
!!!                             |         |
!!!                          source      sink
!!!
!!!
!!!      IF 'key_trc_diabio' key is activated, the biogeochemical
!!!    trends for passive tracers are saved for futher diagnostics.
!!!
!!!      multitasked on vertical slab (jj-loop)
!!!
!!!   MODIFICATIONS:
!!!   --------------

      USE myalloc
      USE BIO_mem
      USE OPT_mem, ONLY: PAR, RMU
      USE BC_mem
      USE mpi

! ----------------------------------------------------------------------
!  BEGIN BC_REFACTORING SECTION
!  ---------------------------------------------------------------------

      use bc_set_mod

! ----------------------------------------------------------------------
!  END BC_REFACTORING SECTION
!  ---------------------------------------------------------------------
      
      IMPLICIT NONE


!!!----------------------------------------------------------------------
!!! local declarations
!!! ==================

      double precision,dimension(jptra,jpk*jpj*jpi) :: b
      double precision,dimension(jpk*jpj*jpi,jptra) :: a
      double precision,dimension(4,jpk*jpj*jpi) :: c
      double precision,dimension(jptra_dia,jpk*jpj*jpi) :: d
      double precision,dimension(jpk*jpj*jpi,16) :: er
      double precision,dimension(jptra_dia_2d) :: d2


      integer :: jk,jj,ji,jb,jn
      integer :: jtr,jtrmax,tra_idx
      integer :: bottom
      double precision :: correct_fact


!!!----------------------------------------------------------------------
!!! statement functions
!!! ===================


!   | --------------|
!   | BFM MODEL CALL|
!   | --------------|

        BIOparttime = MPI_WTIME()

          surf_mask(:) = 0.
          surf_mask(1) = 1.
! -------------------------------------------------

          ! tra_idx = tra_matrix_gib(1)
          jtrmax=jptra

! ---------------- Fuori dai punti BFM

      ogstm_sediPI=0.
      tra_DIA    = 0.
      tra_DIA_2d = 0. ! da sistemare


!    Initialization
      a        = 1.0
      er       = 1.0
      er(:,10) = 8.1


#ifdef gdept1d
! er(:,11) is calculated outside the loop on ji,jj
      do jk=1, jpk
         correct_fact= 1.0D0

         if ( (gdept(jk) .GT. 1000.0D0 ) .AND.  (gdept(jk) .LT. 2000.0D0 )) then
             correct_fact= 0.25D0
         endif

         if (gdept(jk) .GE. 2000.0D0 ) then
             correct_fact= 0.0D0
         endif

         er(jk,11) = correct_fact * ( gdept(jpk)-gdept(jk) ) /gdept(jpk)
      enddo
#endif
       DO jtr=1, jtrmax
       counter=1
       DO ji=1,jpi
       DO jj=1,jpj
       if (bfmmask(1,jj,ji) == 0) then
         counter = counter + jpk
         cycle
       endif
       DO jk=1,jpk
        a(counter, jtr) = trn(jk,jj,ji,jtr) ! current biogeochemical concentrations
        counter=counter+1
       ENDDO
       ENDDO
       ENDDO
       ENDDO

      counter=1
      DO ji=1,jpi
          DO jj=1,jpj
             if (bfmmask(1,jj,ji) == 0) then
                         counter = counter +1
             else

!                bottom = mbathy(jj,ji)
                 DO jtr=1, jtrmax
                  DO jk=1, jpk

                      a(counter, jtr) = trn(jk,jj,ji,jtr) ! current biogeochemical concentrations


! Environmental regulating factors (er,:)

                          er(counter,1)  = tn (jk,jj,ji)! Temperature (Celsius)
                          er(counter,2)  = sn (jk,jj,ji)  ! Salinity PSU
                          er(counter,3)  = rho(jk,jj,ji)        ! Density Kg/m3
                          er(counter ,4)  = ice                  ! from 0 to 1 adimensional
                          if (jk==1) then
                               er(counter ,5)  = ogstm_co2(jj,ji)     ! CO2 Mixing Ratios (ppm)  390
                          else
                               er(counter ,5)  = 0.0D0
                          endif
                          er(counter,6)  = PAR(jk,jj,ji,1) ! PAR for diatoms
                          er(counter,7)  = PAR(jk,jj,ji,2) ! PAR for flagellates
                          er(counter,8)  = PAR(jk,jj,ji,3) ! PAR for pico phytoplankton
                          er(counter,9)  = PAR(jk,jj,ji,4) ! PAR for dinoflagellates
                          er(counter,10) = PAR(jk,jj,ji,5) ! total PAR for CDOM
                          if (jk==1) then
                             er(counter       ,11)  = DAY_LENGTH(jj,ji)    ! fotoperiod expressed in hours
                          else
                             er(counter       ,11)  = 0.0D0
                          endif
                          er(counter,12)  = e3t(jk,jj,ji)        ! depth in meters of the given cell
                          if (jk==1) then
                              er(counter       ,13)  = vatm(jj,ji)                ! wind speed (m/s)
                          else
                              er(counter       ,13)  = 0.0D0
                          endif
                          er(counter,14) = ogstm_PH(jk,jj,ji)   ! 8.1
                          if (jk==1) then
                              er(counter       ,15)  = RMU(jj,ji)                ! avg. cosine direct
                          else
                              er(counter       ,15)  = 0.0D0
                          endif

#ifndef gdept1d
!                        do jk=1,bottom
                             correct_fact= 1.0D0
                             if ( (gdept(jk,jj,ji) .GT. 1000.0D0 ) .AND.  (gdept(jk,jj,ji) .LT. 2000.0D0)) then
                                 correct_fact= 0.25D0
                             endif

                             if (gdept(jk,jj,ji) .GE. 2000.0D0 ) then
                                 correct_fact= 0.0D0
                             endif

                             er(counter,16) = correct_fact * ( gdept(jpk,jj,ji)-gdept(jk,jj,ji) ) /gdept(jpk,jj,ji)
!                        enddo
#endif
                         counter = counter +1

                          END DO
                          END DO
                          END DO
                          END DO

                          call BFM1D_Input_EcologyDynamics(bottom,a,jtrmax,er)

                         call BFM1D_reset()

                         call EcologyDynamics()

                         call BFM1D_Output_EcologyDynamics(b, c, d, d2)

          counter =1               
          DO ji=1,jpi
             DO jj=1,jpj
                   if (bfmmask(1,jj,ji) == 0) then
                         counter = counter +1
                    else

!                bottom = mbathy(jj,ji)
                    DO jk=1, jpk

                          DO jtr=1, jtrmax
                             tra(jk,jj,ji,jtr) = tra(jk,jj,ji,jtr) + b(jtr,counter) ! trend
                             tra_DIA(jtr, jk ,jj,ji) = d(jtr,counter) ! diagnostic
                             ogstm_PH(jk,jj,ji) = d(pppH,counter) ! Follows solver guess, put 8.0 if pppH is not defined
                          END DO

                          DO jtr=1,4
                             ogstm_sediPI(jk,jj,ji,jtr) = c(jtr,counter)      ! BFM output of sedimentation speed (m/d)
                          END DO



                         if (jk==1) then
                            tra_DIA_2d(:,jj,ji) = d2(:) ! diagnostic
                         endif


                         counter = counter +1

                          END DO
                       END DO
                    END DO
                END DO

       DO jtr=1, jtrmax
       counter =1
       DO ji=1,jpi
       DO jj=1,jpj
       if (bfmmask(1,jj,ji) == 0) then
         counter = counter + jpk
         cycle
       endif
       DO jk=1,jpk
          tra(jk,jj,ji,jtr) = tra(jk,jj,ji,jtr) + b(jtr,counter) ! trend
          tra_DIA(jtr, jk ,jj,ji) = d(jtr,counter) ! diagnostic
          ogstm_PH(jk,jj,ji) = d(pppH,counter) ! Follows solver guess, put 8.0 if pppH is not defined
          counter=counter+1
       ENDDO
       ENDDO
       ENDDO
       ENDDO
       DO jtr=1, 4
       counter =1
       DO ji=1,jpi
       DO jj=1,jpj
       if (bfmmask(1,jj,ji) == 0) then
         counter = counter + jpk
         cycle
       endif
       DO jk=1,jpk
       tra(jk,jj,ji,jtr) = tra(jk,jj,ji,jtr) + b(jtr,counter) ! trend
        counter=counter+1
       ENDDO
       ENDDO
       ENDDO
       ENDDO

! ----------------------------------------------------------------------
!  BEGIN BC_REFACTORING SECTION
!  ---------------------------------------------------------------------

      call boundaries%fix_diagnostic_vars(tra_DIA, tra_DIA_2d)

! ----------------------------------------------------------------------
!  END BC_REFACTORING SECTION
!  ---------------------------------------------------------------------


                BIOparttime =  MPI_WTIME() -BIOparttime
                BIOtottime  = BIOtottime  + BIOparttime
               
      END SUBROUTINE trcbio
