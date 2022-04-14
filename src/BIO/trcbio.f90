
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

      integer                    :: counter,counterS,counterB
      integer,dimension(jpj*jpi) :: SRFidx,BOTidx
      double precision,dimension(jpk*jpj*jpi,jptra) :: a,b
      double precision,dimension(jpk*jpj*jpi,4) :: c
      double precision,dimension(jptra_dia,jpk*jpj*jpi) :: d
      double precision,dimension(jpk*jpj*jpi,11) :: er
      double precision,dimension(jptra_dia_2d,jpj*jpi) :: d2


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

       DO jtr=1, jtrmax
       counter=1
       DO ji=1,jpi
       DO jj=1,jpj
       if (bfmmask(1,jj,ji) == 0) then
         counter = counter + jpk
         cycle
       endif
       bottom = mbathy(jj,ji)
       DO jk=1,bottom
        a(counter, jtr) = trn(jk,jj,ji,jtr) ! current biogeochemical concentrations
        counter=counter+1
       ENDDO
       DO jk=bottom+1,jpk
        counter=counter+1
       ENDDO
       ENDDO
       ENDDO
       ENDDO

! Environmental regulating factors (er,:)
       counter=1
       counterS=1
       counterB=1
       DO ji=1,jpi
       DO jj=1,jpj
       if (bfmmask(1,jj,ji) == 0) then
         SRFidx(counterS)=counter
         counterS = counterS + 1
         counter = counter + jpk
         BOTidx(counterB)=counter-1
         counterB = counterB + 1
         cycle
       endif
       bottom = mbathy(jj,ji)
       SRFidx(counterS)=counter
       counterS = counterS + 1
       DO jk=1,bottom
! Environmental regulating factors (er,:)
               er(counter,1)  = tn (jk,jj,ji)  ! Temperature (Celsius)
               er(counter,2)  = sn (jk,jj,ji)  ! Salinity PSU
               er(counter,3)  = rho(jk,jj,ji)  ! Density Kg/m3
               er(counter,4)  = ice            ! from 0 to 1 adimensional

               if (jk==1) then
                    er(counter ,5)  = ogstm_co2(jj,ji)     ! CO2 Mixing Ratios (ppm)  390
               else
                    er(counter ,5)  = huge(er(1,1))
               endif

               er(counter,6) = instant_par(COMMON_DATEstring,xpar(jk,jj,ji))  ! PAR umoles/m2/s | Watt to umoles photons W2E=1./0.217

               if (jk==1) then
                    er(counter   ,7)  = DAY_LENGTH(jj,ji)    ! fotoperiod expressed in hours
               else
                    er(counter   ,7)  = huge(er(1,1))
               endif

               er(counter,8)  = e3t(jk,jj,ji)        ! depth in meters of the given cell

               if (jk==1) then
                    er(counter,9)  = vatm(jj,ji)                ! wind speed (m/s)
               else
                    er(counter,9)  = huge(er(1,1))
               endif

               er(counter,10) = ogstm_PH(jk,jj,ji)   ! 8.1

               correct_fact= 1.0D0
               if ( (gdept(jk,jj,ji) .GT. 1000.0D0 ) .AND.  (gdept(jk,jj,ji) .LT. 2000.0D0)) then
                    correct_fact= 0.25D0
               endif

               if (gdept(jk,jj,ji) .GE. 2000.0D0 ) then
                    correct_fact= 0.0D0
               endif

               er(counter,11) = correct_fact * ( gdept(jpk,jj,ji)-gdept(jk,jj,ji) ) /gdept(jpk,jj,ji)

               counter=counter+1

       ENDDO
       BOTidx(counterB)=counter-1
       counterB = counterB + 1
       DO jk=bottom+1,jpk
        counter=counter+1
       ENDDO
       ENDDO
       ENDDO

!      call BFM1D_Input_EcologyDynamics(bottom,a,jtrmax,er)

       call BFM1D_Input_EcologyDynamics(SRFidx,BOTidx,a,jtrmax,er)

       call BFM1D_reset()

       call EcologyDynamics()

       call BFM1D_Output_EcologyDynamics(b, c, d, d2)

       DO jtr=1, jtrmax
       counter =1
       DO ji=1,jpi
       DO jj=1,jpj
       if (bfmmask(1,jj,ji) == 0) then
            counter = counter + jpk
            cycle
       endif

       DO jk=1, jpk
          tra(jk,jj,ji,jtr) = tra(jk,jj,ji,jtr) + b(counter,jtr) * bfmmask(jk,jj,ji) ! trend
          counter = counter + 1
       ENDDO

       ENDDO
       ENDDO
       ENDDO

       counter=1
       DO ji=1,jpi
       DO jj=1,jpj
       if (bfmmask(1,jj,ji) == 0) then
            counter = counter + jpk
            cycle
       endif

       DO jk=1, jpk

       DO jtr=1, jtrmax
          tra_DIA(jtr, jk ,jj,ji) = d(jtr,counter) * bfmmask(jk,jj,ji) ! diagnostic
       ENDDO
          ogstm_PH(jk,jj,ji) = d(pppH,counter) ! Follows solver guess, put 8.0 if pppH is not defined
       DO jtr=1, 4
          ogstm_sediPI(jk,jj,ji,jtr) = c(counter,jtr) * bfmmask(jk,jj,ji)     ! BFM output of sedimentation speed (m/d)
       ENDDO

          counter = counter + 1
       ENDDO

       ENDDO
       ENDDO

       counter=1
       DO ji=1,jpi
       DO jj=1,jpj
           tra_DIA_2d(:,jj,ji) = d2(:,counter) * bfmmask(1,jj,ji)! diagnostic
           counter = counter +1
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
