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

  integer :: jk, jj, ji, jn, jlinear2d, jlinear3d, bottom, queue
  double precision :: correct_fact, gdept_local, gdeptmax_local
  
  integer :: year, month, day
  double precision :: sec

  BIOparttime = MPI_WTIME()
  
  ! Initialization

  queue=1

  !$acc kernels default(present) async(queue)
  D3STATE = 1.0
  er = 1.0
  er(:,10) = 8.1
  tra_DIA = 0.
  tra_DIA_2d = 0. ! da sistemare (?)

  ! ogstm_sediPI appear to be unused
  ogstm_sediPI = 0.
  !$acc end kernels

  ! NOTE: this kernel *needs* to be executed synchronously as we need the
  ! reduced `bottom` scalar value on host before launching the next kernel.
  bottom=0
  !$acc parallel loop gang vector reduction(max:bottom) collapse(2) default(present)
  do ji = 1, jpi
     do jj = 1, jpj
        bottom = max(bottom,mbathy(jj,ji))
     end do
  end do
  !$acc end parallel loop

  ! Set D3STATE (pass state to BFM)
  !$acc parallel loop gang vector default(present) collapse(4) async(queue)
  do jn = 1, jptra
    do ji = 1, jpi
      do jj = 1, jpj
         do jk = 1, bottom
            if (.not. bfmmask(1, jj, ji) == 0) then
               if (jk <= mbathy(jj, ji)) then
                  jlinear3d = jk + (jj - 1) * jpk + (ji - 1) * jpk * jpj
                  D3STATE(jlinear3d, jn) = trn(jk, jj, ji, jn)
               endif
            endif
          end do
      end do
    end do
  end do
  !$acc end parallel loop

  call read_date_string(COMMON_DATEstring, year, month, day, sec)

  ! Set er
  !$acc parallel loop gang vector default(present) collapse(3) async(queue)
  do ji = 1, jpi
     do jj = 1, jpj
        do jk = 1, bottom
           if (.not. bfmmask(1, jj, ji) == 0) then
              if (jk <= mbathy(jj, ji)) then
                 jlinear3d = jk + (jj - 1) * jpk + (ji - 1) * jpk * jpj
                 if (jk .eq. 1) then
                    er(jlinear3d, 4) = ice ! from 0 to 1 adimensional
                    er(jlinear3d, 5) = ogstm_co2(jj, ji) ! CO2 Mixing Ratios (ppm)  390
                    er(jlinear3d, 7) = DAY_LENGTH(jj, ji) ! fotoperiod expressed in hours
                    er(jlinear3d, 9) = vatm(jj, ji) ! wind speed (m/s)
                 end if
                 er(jlinear3d, 1) = tn(jk, jj, ji) ! Temperature (Celsius)
                 er(jlinear3d, 2) = sn(jk, jj, ji) ! Salinity PSU
                 er(jlinear3d, 3) = rho(jk, jj, ji) ! Density Kg/m3
                 er(jlinear3d, 6) = instant_par_from_sec(sec, xpar(jk, jj, ji)) ! PAR umoles/m2/s | Watt to umoles photons W2E=1./0.217
                 er(jlinear3d, 8) = e3t(jk, jj, ji) ! depth in meters of the given cell
                 er(jlinear3d, 10) = ogstm_PH(jk, jj, ji) ! 8.1
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
                 er(jlinear3d, 11) = correct_fact * (gdeptmax_local - gdept_local) / gdept_local
              endif
           end if
        end do
     end do
  end do
  !$acc end parallel loop

  !$acc wait(queue)
  !$acc update host(D3STATE)

  call BFM1D_Input_EcologyDynamics(mbathy, er) ! here mbathy was bottom
  call BFM1D_reset()
  call EcologyDynamics()
  call BFM1D_Output_EcologyDynamics(sediPPY, local_D3DIAGNOS, local_D2DIAGNOS)

  !$acc update device(D3SOURCE,sediPPY,local_D3DIAGNOS,local_D2DIAGNOS)

  !$acc parallel loop gang vector collapse(4) default(present) async(queue)
  do jn = 1, max(4, jptra, jptra_dia)
    do ji = 1, jpi
      do jj = 1, jpj
        do jk = 1, bottom
          if (.not. bfmmask(1, jj, ji) == 0) then
            if (jk <= mbathy(jj, ji)) then
               jlinear3d = jk + (jj - 1) * jpk + (ji - 1) * jpk * jpj
               if (jn .le. jptra) then
                  tra(jk, jj, ji, jn) = tra(jk, jj, ji, jn) + D3SOURCE(jlinear3d, jn) ! trend
               endif
               if (jn .le. jptra_dia) then
                  tra_DIA(jk, jj, ji, jn) = local_D3DIAGNOS(jlinear3d, jn)
               endif
               if (jn .le. 4) then
                  ogstm_sediPI(jk, jj, ji, jn) = sediPPY(jlinear3d, jn)   ! BFM output of sedimentation speed (m/d)
               endif
            end if
          end if
        end do
      end do
    end do
  end do
  !$acc end parallel loop

  !$acc parallel loop gang vector collapse(3) default(present) async(queue)
  do ji = 1, jpi
    do jj = 1, jpj
      do jk = 1, bottom
         if (.not. bfmmask(1, jj, ji) == 0) then
            if (jk <= mbathy(jj, ji)) then
               jlinear2d = jj + (ji - 1) * jpj
               tra_DIA_2d(:, jj, ji) = local_D2DIAGNOS(jlinear2d, :)
               jlinear3d = jk + (jj - 1) * jpk + (ji - 1) * jpk * jpj
               ogstm_PH(jk, jj, ji) = local_D3DIAGNOS(jlinear3d, pppH) ! Follows solver guess, put 8.0 if pppH is not defined
            end if
         end if
      end do
    end do
  end do
  !$acc end parallel loop

  !---------------------------------------------------------------------
  ! BEGIN BC_REFACTORING SECTION
  !---------------------------------------------------------------------
  ! XXX: when should we care about this ?
  call boundaries%fix_diagnostic_vars()
  !----------------------------------------------------------------------
  ! END BC_REFACTORING SECTION
  !---------------------------------------------------------------------

  !$acc wait(queue)

  BIOparttime = MPI_WTIME() - BIOparttime
  BIOtottime = BIOtottime + BIOparttime
END SUBROUTINE trcbio
