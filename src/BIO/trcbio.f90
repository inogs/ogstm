
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
      ! epascolo USE myalloc_mpp
      USE BIO_mem
      USE BC_mem
      use global_mem, ONLY:RLEN 
      USE mpi
      
      IMPLICIT NONE


!!!----------------------------------------------------------------------
!!! local declarations
!!! ==================
      logical :: sur,bot
      double precision,dimension(:,:,:,:),allocatable :: A
      double precision,dimension(jptra) :: b
      double precision,dimension(4) :: c
      double precision,dimension(jptra_dia) :: d
      double precision,dimension(10) :: er
      double precision,dimension(jptra_dia_2d) :: d2

      integer :: jk,jj,ji,jb,jn
      integer :: jtr,jtrmax,tra_idx


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

          tra_idx = tra_matrix_gib(1)
          jtrmax=jptra

! ---------------- Fuori dai punti BFM

      sediPI     = 0.
      tra_DIA    = 0.
      tra_DIA_2d = 0.

allocate(A(jptra,jpk,jpj,jpi))

 DO jn=1,jptra
  DO ji=1,jpi
   DO jj=1,jpj
    DO jk=1,jpk
      A(jn,jk,jj,ji) =  trn(jk,jj,ji,jn)
    END DO
   END DO
  END DO
 END DO



! $omp   parallel do default(none)  private(jb,jk,jj,ji,mytid,sur,bot,jtr,a,b,c,d,d2,er)
! $omp&      shared(NBFMPOINTS, BFMpoints,tra_idx,tra_matrix_gib,
! $omp&               restotr,jtrmax,trn,tn,sn,xpar,e3t,vatm,surf_mask,DAY_LENGTH,
! $omp&             sediPI,PH,tra_DIA,tra_DIA_2d,tra,rho,ice,co2,idxt2glo)

      MAIN_LOOP: DO  jb = 1, NBFMPOINTS

                 !IF( mytid + jb <= NBFMPOINTS ) THEN


                 ji = BFMpoints(3, jb)
                 jj = BFMpoints(2, jb)
                 jk = BFMpoints(1, jb)

                        IF(jk .eq. 1) then
                              SRFindices = .TRUE.
                        ELSE
                              SRFindices = .FALSE.
                        ENDIF

                        bot = .FALSE.
          
                       D3STATE(:) = A(:,jk,jj,ji)

                        O2o = A(ppO2o,jk,jj,ji)
                        N1p = A(ppN1p,jk,jj,ji)
                        N3n = A(ppN3n,jk,jj,ji)
                        N4n = A(ppN4n,jk,jj,ji)
                        O4n = A(ppO4n,jk,jj,ji)
                        N5s = A(ppN5s,jk,jj,ji)
                        N6r = A(ppN6r,jk,jj,ji)
                        B1c = A(ppB1c,jk,jj,ji)
                        B1n = A(ppB1n,jk,jj,ji)
                        B1p = A(ppB1p,jk,jj,ji)
                        P1c = A(ppP1c,jk,jj,ji)
                        P2c = A(ppP2c,jk,jj,ji)
                        P3c = A(ppP3c,jk,jj,ji)
                        P4c = A(ppP4c,jk,jj,ji)
                        P1n = A(ppP1n,jk,jj,ji)
                        P2n = A(ppP2n,jk,jj,ji)
                        P3n = A(ppP3n,jk,jj,ji)
                        P4n = A(ppP4n,jk,jj,ji)
                        P1p = A(ppP1p,jk,jj,ji)
                        P2p = A(ppP2p,jk,jj,ji)
                        P3p = A(ppP3p,jk,jj,ji)
                        P4p = A(ppP4p,jk,jj,ji)
                        P1s = A(ppP1s,jk,jj,ji)
                        P1i = A(ppP1l,jk,jj,ji)
                        P2i = A(ppP2l,jk,jj,ji)
                        P3i = A(ppP3l,jk,jj,ji)
                        P4i = A(ppP4l,jk,jj,ji)
                        Z3c = A(ppZ3c,jk,jj,ji)
                        Z4c = A(ppZ4c,jk,jj,ji)
                        Z3n = A(ppZ3n,jk,jj,ji)
                        Z4n = A(ppZ4n,jk,jj,ji)
                        Z3p = A(ppZ3p,jk,jj,ji)
                        Z4p = A(ppZ4p,jk,jj,ji)
                        Z5c = A(ppZ5c,jk,jj,ji)
                        Z6c = A(ppZ6c,jk,jj,ji)
                        Z5n = A(ppZ5n,jk,jj,ji)
                        Z6n = A(ppZ6n,jk,jj,ji)
                        Z5p = A(ppZ5p,jk,jj,ji)
                        Z6p = A(ppZ6p,jk,jj,ji)
                        R1c = A(ppR1c,jk,jj,ji)
                        R1n = A(ppR1n,jk,jj,ji)
                        R1p = A(ppR1p,jk,jj,ji)
                        R1s = A(ppR1s,jk,jj,ji)
                        R2c = A(ppR2c,jk,jj,ji)
                        R6c = A(ppR6c,jk,jj,ji)
                        R6n = A(ppR6n,jk,jj,ji)
                        R6p = A(ppR6p,jk,jj,ji)
                        R6s = A(ppR6s,jk,jj,ji)
                        R7c = A(ppR7c,jk,jj,ji)
#ifdef  INCLUDE_PELCO2
                        O3c = A(ppO3c,jk,jj,ji)
                        O3h = A(ppO3h,jk,jj,ji)
#endif

                        !   DO jtr=1, jtrmax
                        !      a(jtr) = trn(jk,jj,ji,jtr) ! current biogeochemical concentrations
                        !       WRITE(*,200) ,'I',jk,jj,ji,jtr,trn(jk,jj,ji,jtr)   
                        !   END DO
! Environmental regulating factors (er)

                          ETW  = tn (jk,jj,ji)        ! Temperature (Celsius)
                          ESW  = sn (jk,jj,ji)        ! Salinity PSU
                          ERHO  = rho(jk,jj,ji)        ! Density Kg/m3
                          EICE  = ice                  ! from 0 to 1 adimensional
                          EPCO2air  = ogstm_co2(jj,ji)           ! CO2 Mixing Ratios (ppm)  390
                          BFM0D_Irr  = xpar(jk,jj,ji)*86400.0_RLEN      ! PAR umoles/m2/s | Watt to umoles photons W2E=1./0.217
                          SUNQ  = DAY_LENGTH(jj,ji)    ! fotoperiod expressed in hours
                          DEPTH  = e3t(jk,jj,ji)        ! depth in meters of the given cell
                          EWIND  = vatm(jj,ji) * surf_mask(jk) ! wind speed (m/s)
                          ph = ogstm_PH(jk,jj,ji)         ! PH
                        !   WRITE(*,201),'ERR',jk,jj,ji,er(:)
                          !call BFM0D_Input_EcologyDynamics(sur,bot,A(:,jk,jj,ji),jtrmax,er)

                          call BFM0D_reset()

                         call EcologyDynamics()

                        !   if (SRFindices) then
                        !      call BFM0D_Output_EcologyDynamics_surf(b, c, d ,d2)
                        !    else
                        !       call BFM0D_Output_EcologyDynamics(b, c, d)
                        !   endif

                          tra(jk,jj,ji,:) =tra(jk,jj,ji,:) + D3SOURCE(:)/86400.0_RLEN

                        !   DO jtr=1, jtrmax
                        !      tra(jk,jj,ji,jtr) =tra(jk,jj,ji,jtr) +b(jtr) ! trend
                        ! !      WRITE(*,200),'OB',jk,jj,ji,jtr,b(jtr)
                        ! !      WRITE(*,200),'TRA',jk,jj,ji,jtr,tra(jk,jj,ji,jtr)
                        !   END DO

                          !! sedimenti
                          ogstm_sediPI(jk,jj,ji,1) = sediPI( iiP1 ) ! BFM output of sedimentation speed (m/d)
                          ogstm_sediPI(jk,jj,ji,2) = sediPI( iiP2 ) ! BFM output of sedimentation speed (m/d)
                          ogstm_sediPI(jk,jj,ji,3) = sediPI( iiP3 ) ! BFM output of sedimentation speed (m/d)
                          ogstm_sediPI(jk,jj,ji,4) = sediPI( iiP4 ) ! BFM output of sedimentation speed (m/d)

                        
                        tra_DIA(jk,jj,ji,1)  = BFM0D_ppg(1) + BFM0D_ppg(2) + BFM0D_ppg(3) + BFM0D_ppg(4) ! gross primary production:mg C/m3/d
                        tra_DIA(jk,jj,ji,2)  = BFM0D_ppn(1) + BFM0D_ppn(2) + BFM0D_ppn(3) + BFM0D_ppn(4) ! net primary production:mg C/m3/d
                        tra_DIA(jk,jj,ji,3)  = BFM0D_ppb ! bacterial production:mg C/m3/d

                        !CO2 module state variables

#ifdef INCLUDE_PELCO2
                        tra_DIA(jk,jj,ji,4)  = DIC
                        tra_DIA(jk,jj,ji,5)  = CO2
                        tra_DIA(jk,jj,ji,6)  = HCO3
                        tra_DIA(jk,jj,ji,7)  = CO3
                        tra_DIA(jk,jj,ji,8)  = Ac
                        tra_DIA(jk,jj,ji,9)  = pH
                        tra_DIA(jk,jj,ji,10) = pCO2
                        tra_DIA(jk,jj,ji,11) = - D3SOURCE(ppN3n) - D3SOURCE(ppN1p)! impact on alcalinity of NO3 and PO4 uptake (F02)
                        tra_DIA(jk,jj,ji,12) =   D3SOURCE(ppN4n)                  ! impact of ammonia uptake                   (F03)
#endif
                        tra_DIA(jk,jj,ji,13) = BFM0D_phr(1) ! respiration diatoms             (F04)
                        tra_DIA(jk,jj,ji,14) = BFM0D_phr(2) ! respiration flagellates         (F05)
                        tra_DIA(jk,jj,ji,15) = BFM0D_phr(3) ! respiration pico                (F06)
                        tra_DIA(jk,jj,ji,16) = BFM0D_phr(4) ! respiration dinaflagellates     (F07)
                        tra_DIA(jk,jj,ji,17) = BFM0D_bar    ! respiration bacteria            (F08)
                        tra_DIA(jk,jj,ji,18) = BFM0D_uzr(1) ! respiration microzooplancton    (F09)
                        tra_DIA(jk,jj,ji,19) = BFM0D_uzr(2) ! respiration nanoflagellates     (F10)
                        tra_DIA(jk,jj,ji,20) = BFM0D_mzr(1) ! respiration meso carnivorous    (F11)
                        tra_DIA(jk,jj,ji,21) = BFM0D_mzr(2) ! respiration meso omnivorous     (F12)
                        if (SRFindices) then
                        tra_DIA_2d(jj,ji,1) =CO2airflux
                        endif
                        !   DO jtr=1,jptra_dia
                        !      tra_DIA(jk,jj,ji,jtr) = d(jtr) ! diagnostic
                        ! !      WRITE(*,200),'IO3',jk,jj,ji,jtr,tra_DIA(jk,jj,ji,jtr)
                        !   END DO
                          
                        !   if (SRFindices) then
                        !       DO jtr=1,jptra_dia_2d
                        !          tra_DIA_2d(jj,ji,jtr) = d2(jtr) ! diagnostic
                        !       !    WRITE(*,202),'IO2',jj,ji,jtr,tra_DIA_2d(jj,ji,jtr)
                        !       END DO
                        !   endif

                          ogstm_PH(jk,jj,ji)=pH ! Follows solver guess, put 8.0 if pppH is not defined


             !ENDIF

                END DO MAIN_LOOP
! 200    FORMAT(' ',A3,I4,I4,I4,I4,D30.23)
! 201    FORMAT(' ',A3,I4,I4,I4,10D30.23)
! 202    FORMAT(' ',A3,I4,I4,I4,10D30.23)
                

                  deallocate(A)                          
! $omp end parallel do

                BIOparttime =  MPI_WTIME() -BIOparttime
                BIOtottime  = BIOtottime  + BIOparttime
               
      END SUBROUTINE trcbio
