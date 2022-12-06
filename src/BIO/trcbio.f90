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
SUBROUTINE trcbio
  use myalloc
  use BIO_mem
  use BC_mem
  use mpi
  use mem, only: D3STATE, D3SOURCE

  !----------------------------------------------------------------------
  ! BEGIN BC_REFACTORING SECTION
  ! ---------------------------------------------------------------------
  use bc_set_mod
  !----------------------------------------------------------------------
  ! END BC_REFACTORING SECTION
  ! ---------------------------------------------------------------------

  IMPLICIT NONE

  double precision, dimension(jpi * jpj * jpk, 4) :: sediPPY
  double precision, dimension(jpi * jpj * jpk, jptra_dia) :: local_D3DIAGNOS
  double precision, dimension(jptra_dia_2d) :: local_D2DIAGNOS
  double precision, dimension(jpi * jpj * jpk, 11) :: er

  integer :: jk, jj, ji, jn, jlinear, bottom
  double precision :: correct_fact, gdept_local, gdeptmax_local

  BIOparttime = MPI_WTIME()
  
  ! Initialization
  D3STATE = 1.0
  er = 1.0
  er(:,10) = 8.1
  tra_DIA = 0.
  tra_DIA_2d = 0. ! da sistemare (?)

  ! ogstm_sediPI appear to be unused
  ogstm_sediPI = 0.

  ! Set D3STATE (pass state to BFM)
  do jn = 1, jptra
    do ji = 1, jpi
      do jj = 1, jpj
        if (bfmmask(1, jj, ji) == 0) then
          cycle
        else
          bottom = mbathy(jj, ji)
          do jk = 1, bottom
            jlinear = jk + (jj - 1) * jpk + (ji - 1) * jpk * jpj
            D3STATE(jlinear, jn) = trn(jk, jj, ji, jn)
          end do
        end if
      end do
    end do
  end do

  ! Set er 
  do ji = 1, jpi
    do jj = 1, jpj
      if (bfmmask(1, jj, ji) == 0) then
        cycle
      else
        bottom = mbathy(jj, ji)
        do jk = 1, bottom
          jlinear = jk + (jj - 1) * jpk + (ji - 1) * jpk * jpj
          if (jk .eq. 1) then
            er(jlinear, 4) = ice ! from 0 to 1 adimensional
            er(jlinear, 5) = ogstm_co2(jj, ji) ! CO2 Mixing Ratios (ppm)  390
            er(jlinear, 7) = DAY_LENGTH(jj, ji) ! fotoperiod expressed in hours
            er(jlinear, 9) = vatm(jj, ji) ! wind speed (m/s)
          end if
          er(jlinear, 1) = tn(jk, jj, ji) ! Temperature (Celsius)
          er(jlinear, 2) = sn(jk, jj, ji) ! Salinity PSU
          er(jlinear, 3) = rho(jk, jj, ji) ! Density Kg/m3
          er(jlinear, 6) = instant_par(COMMON_DATEstring, xpar(jk, jj, ji)) ! PAR umoles/m2/s | Watt to umoles photons W2E=1./0.217
          er(jlinear, 8) = e3t(jk, jj, ji) ! depth in meters of the given cell
          er(jlinear, 10) = ogstm_PH(jk, jj, ji) ! 8.1
          #ifdef gdept1d
          gdept_local = gdept(jk)
          gdeptmax_local = gdept(jpk)
          #else
          gdept_local = gdept(jk, jj, ji)
          gdeptmax_local = gdept(jpk, jj, ji)
          #endif
          if (gdept_local .lt. 1000.0D0) then
            correct_fact = 1.0D0
          else if (gdept_local .lt. 2000.0D0) then
            correct_fact = 0.25D0
          else
            correct_fact = 0.0D0
          end if
          er(jlinear, 11) = correct_fact * (gdeptmax_local - gdept_local) / gdept_local
        end do
      end if
    end do
  end do
   
  write (*, *) "TRA surf: ", tra(1, :, :, 49)
  write (*, *) "TRA min, argmin: ", minval(tra(:, :, :, 49)), minloc(tra(:, :, :, 49))
  write (*, *) "d3source min, argmin: ", minval(d3source(:, 49)), minloc(d3source(:, 49))
  write (*, *) "d3state min, argmin: ", minval(d3state(:, 49)), minloc(d3state(:, 49))
  call BFM1D_Input_EcologyDynamics(mbathy, er) ! here mbathy was bottom
  call BFM1D_reset()
  call EcologyDynamics()
  call BFM1D_Output_EcologyDynamics(sediPPY, local_D3DIAGNOS, local_D2DIAGNOS)
  write (*, *) "d3source min, argmin: ", minval(d3source(:, 49)), minloc(d3source(:, 49))
  write (*, *) "d3state min, argmin: ", minval(d3state(:, 49)), minloc(d3state(:, 49))

  do jn = 1, jptra ! max(4, jptra, jptra_dia)
    do ji = 1, jpi
      do jj = 1, jpj
        if (bfmmask(1, jj, ji) == 0) then
          cycle
        else
          bottom = mbathy(jj, ji)
          do jk = 1, bottom
            jlinear = jk + (jj - 1) * jpk + (ji - 1) * jpk * jpj
            !if (jn .le. jptra) then
              tra(jk, jj, ji, jn) = tra(jk, jj, ji, jn) + D3SOURCE(jlinear, jn) ! trend
            !end if
            !if (jn .le. jptra_dia) then 
            !  tra_DIA(jk, jj, ji, jn) = local_D3DIAGNOS(jlinear, jn)
            !end if
            !if (jn .le. 4) then
            !  ogstm_sediPI(jk, jj, ji, jn) = sediPPY(jlinear, jn)   ! BFM output of sedimentation speed (m/d)
            !end if
          end do
        end if
      end do
    end do
  end do
  
  write (*, *) "TRA min, argmin: ", minval(tra(:, :, :, 49)), minloc(tra(:, :, :, 49))

  !do ji = 1, jpi
  !  do jj = 1, jpj
  !    if (bfmmask(1, jj, ji) == 0) then
  !      cycle
  !    else
  !      bottom = mbathy(jj, ji)
  !      tra_DIA_2d(:, jj, ji) = local_D2DIAGNOS(:)
  !      do jk = 1, bottom
  !        jlinear = jk + (jj - 1) * jpk + (ji - 1) * jpk * jpj
  !        ogstm_PH(jk, jj, ji) = local_D3DIAGNOS(jlinear, pppH) ! Follows solver guess, put 8.0 if pppH is not defined
  !      end do
  !    end if
  !  end do
  !end do

  !---------------------------------------------------------------------
  ! BEGIN BC_REFACTORING SECTION
  !---------------------------------------------------------------------
  call boundaries%fix_diagnostic_vars()
  !----------------------------------------------------------------------
  ! END BC_REFACTORING SECTION
  !---------------------------------------------------------------------
  BIOparttime = MPI_WTIME() - BIOparttime
  BIOtottime = BIOtottime + BIOparttime
END SUBROUTINE trcbio
