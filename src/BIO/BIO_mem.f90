       MODULE BIO_mem 

       USE modul_param 
       USE myalloc
       USE TIME_MANAGER

#ifdef key_trc_fabm
       USE fabm
       USE fabm_types
       USE OPT_mem
#endif

#ifdef Mem_Monitor
       USE check_mem
       USE iso_c_binding
#endif


       IMPLICIT NONE

       public

      double precision, allocatable :: bfm_trn(:), bfm_tra(:)
!!!$omp  threadprivate(bfm_trn,bfm_tra)
      double precision, allocatable :: surf_mask(:)
      double precision, allocatable :: ogstm_sedipi(:,:,:,:)
      double precision, allocatable :: ogstm_ph(:,:,:) ! GUESS for FOLLOWS algorithm
      double precision, allocatable :: NPPF2(:,:,:)
      double precision, allocatable :: ogstm_co2(:,:), co2_IO(:,:,:)
      double precision:: ice
#ifdef key_trc_fabm
      class (type_fabm_model), pointer :: model_fabm
#endif


!!!----------------------------------------------------------------------
      CONTAINS

#ifdef key_trc_fabm
      subroutine initialize_FABM()
        ! Provide extents of the spatial domain (number of layers nz for a 1D column)
        model_fabm => fabm_create_model()
        call model_fabm%set_domain(jpk,jpj,jpi,1.0d0)
        ! At this point (after the call to fabm_create_model), memory should be
        ! allocated to hold the values of all size(model%interior_state_variables) state variables.
        ! Where this memory resides and how it is laid out is typically host-specific.
        ! Below, we assume all state variable values are combined in an array interior_state with
        ! shape nx, ny, nz, size(model%interior_state_variables).
        jptra=size(model_fabm%interior_state_variables)

! In FABM interior diagnostics includes already fluxes
! and they are counted in jptra_var therefore jptra_flux = 0
! we keep definition of jptra_flux for back compatibility with older BFM
! code.
        jptra_var=size(model_fabm%interior_diagnostic_variables)
        jptra_flux=0                                     
        jptra_dia=jptra_var+jptra_flux
        jptra_dia_2d=size(model_fabm%horizontal_diagnostic_variables)
! optical modules
      if (is_FABM_instance(model_fabm,'light_atm')) then
         oasim_fabm=.TRUE.
         oasim_bioptimod=.FALSE.
         write(numout,*) 'FABM-OASIM optical model will be used in OGSTM'
      end if
      if (is_FABM_instance(model_fabm,'lightspectral')) then
         oasim_fabm=.FALSE.
         oasim_bioptimod=.TRUE.
         write(numout,*) 'BIOPTIMOD OASIM optical model will be used in OGSTM'
      end if
      END subroutine initialize_FABM
#endif

      subroutine myalloc_BIO()
      USE calendar

      INTEGER  :: err
      INTEGER  :: ivar
      INTEGER  :: ji,jj,jn,jl
      INTEGER  :: year, month, day
      double precision :: sec
      double precision  :: aux_mem
      character(LEN=256) :: varname
#ifdef key_trc_fabm
      type (type_fabm_interior_variable_id)   :: interior_id
      type (type_fabm_horizontal_variable_id) :: horizontal_id
      type (type_fabm_scalar_variable_id)     :: scalar_id,id_yearday 
#endif

#ifdef Mem_Monitor
       aux_mem = get_mem(err)
#endif

!!!$omp parallel default(none)
       allocate(bfm_trn(jptra))        
       bfm_trn   = huge(bfm_trn(1))
       allocate(bfm_tra(jptra))        
       bfm_tra   = huge(bfm_tra(1))
!!!$omp end parallel
       allocate(surf_mask(jpk))        
       surf_mask = huge(surf_mask(1))
       allocate(ogstm_co2(jpj,jpi))          
       ogstm_co2       = huge(ogstm_co2(1,1))
       allocate(co2_IO(jpj,jpi,2))    
        
       co2_IO    = huge(co2_IO(1,1,1))

#ifdef key_trc_bfm
       allocate(ogstm_sedipi(jpk,jpj,jpi,4))
#endif
       
#ifdef key_trc_fabm
       allocate(ogstm_sedipi(jpk,jpj,jpi,jptra))
#endif 
       
       ogstm_sedipi    = huge(ogstm_sedipi(1,1,1,1))
       allocate(ogstm_ph(jpk,jpj,jpi))       
       ogstm_ph        = huge(ogstm_ph(1,1,1))
       ogstm_ph=8.0
       allocate(NPPF2(jpk,jpj,jpi))
       NPPF2 = 0 ! nut huge, because it will be assigned only in trcBIO in BFMpoints
                 ! and used in hard_tissue_pump.F also in land points
       ice=0
#ifdef key_trc_fabm
        ! Provide FABM with the vertical indices of the surface and bottom, and the land-sea mask.
        call model_fabm%set_bottom_index(mbathy)  ! NB mbkt extents should match dimension lengths provided to model%set_domain
        call model_fabm%set_mask(tmask,tmask(1,:,:)) ! NB tmask extents should match dimension lengths provided to model%set_domain

        ! At this point (after the call to fabm_create_model), memory should be
        ! allocated to hold the values of all size(model%interior_state_variables) state variables.
        ! Where this memory resides and how it is laid out is typically host-specific.
        ! Below, we assume all state variable values are combined in an array interior_state with
        ! shape nx, ny, nz, size(model%interior_state_variables).

        ! Point FABM to your state variable data
        do ivar = 1, size(model_fabm%interior_state_variables)
           call model_fabm%link_interior_state_data(ivar, trn(:,:,:,ivar))
        end do
        ! Point FABM to environmental data, here shown for temperature
        ! Array temp with extents nx,ny,nz is assumed to be allocated.
        ! Do this for all variables on FABM's standard variable list that the model can provide.
        ! For this list, visit https://fabm.net/standard_variables

        call model_fabm%link_interior_data(fabm_standard_variables%temperature, tn) !  Celsius
        call model_fabm%link_interior_data(fabm_standard_variables%practical_salinity, sn) ! PSU
        call model_fabm%link_interior_data(fabm_standard_variables%density, rho) ! kg m-3
        call model_fabm%link_interior_data(fabm_standard_variables%pressure, gdept) ! dbar
        call model_fabm%link_horizontal_data(fabm_standard_variables%mole_fraction_of_carbon_dioxide_in_air,  ogstm_co2) ! CO2 Mixing Ratios (ppm)  
        call model_fabm%link_interior_data(fabm_standard_variables%depth, gdept ) ! m  
        call model_fabm%link_interior_data(fabm_standard_variables%cell_thickness, e3t ) ! m  
        call model_fabm%link_horizontal_data(fabm_standard_variables%wind_speed, vatm) ! m/s 
        call model_fabm%link_horizontal_data(fabm_standard_variables%longitude, glamt) ! degree_east
        call model_fabm%link_horizontal_data(fabm_standard_variables%latitude, gphit) ! degree_north 

        horizontal_id = model_fabm%get_horizontal_variable_id('cloud_area_fraction')
        call model_fabm%link_horizontal_data(horizontal_id, tcc) ! [0-1] 


        horizontal_id = model_fabm%get_horizontal_variable_id('atmosphere_mass_content_of_cloud_liquid_water')
        call model_fabm%link_horizontal_data(horizontal_id, tclw) ! [0-1] 


        horizontal_id=model_fabm%get_horizontal_variable_id('surface_downwelling_shortwave_flux')
        ! to be provided in case monospectral formulation is used
        surface_downwelling_shortwave_flux=500.0d0
        call model_fabm%link_horizontal_data(horizontal_id, surface_downwelling_shortwave_flux) ! W m^-2

        horizontal_id = model_fabm%get_horizontal_variable_id('atmosphere_mass_content_of_water_vapor')
        atmosphere_mass_content_of_water_vapor=0.1d0
        call model_fabm%link_horizontal_data(horizontal_id, atmosphere_mass_content_of_water_vapor) ! kg m^-2

         horizontal_id = model_fabm%get_horizontal_variable_id('visibility_in_air')
         visibility_in_air = 25000.d0
         call model_fabm%link_horizontal_data(horizontal_id, visibility_in_air) ! m

        horizontal_id = model_fabm%get_horizontal_variable_id('aerosol_air_mass_type')
        aerosol_air_mass_type = 1.d0
        call model_fabm%link_horizontal_data(horizontal_id, aerosol_air_mass_type) ! -

         horizontal_id = model_fabm%get_horizontal_variable_id('surface_specific_humidity')
         surface_specific_humidity = 0.01d0
        call model_fabm%link_horizontal_data(horizontal_id, surface_specific_humidity) ! kg kg^-1

        horizontal_id = model_fabm%get_horizontal_variable_id('surface_temperature')
        call model_fabm%link_horizontal_data(horizontal_id, t2m) ! - degree_Celsius

        horizontal_id = model_fabm%get_horizontal_variable_id('surface_air_pressure')
        call model_fabm%link_horizontal_data(horizontal_id, sp) ! - Pa

       id_yearday = model_fabm%get_scalar_variable_id(fabm_standard_variables%number_of_days_since_start_of_the_year)
       yearday = 1.5d0
       call model_fabm%link_scalar(id_yearday, yearday) ! - days

      if (oasim_bioptimod) then
            do jl=1,nlt
                  
                  varname = 'surf_direct_downward_irradiance_' // lam_strings(jl) // '_nm'
                  horizontal_id = model_fabm%get_horizontal_variable_id(TRIM(varname))
                  call model_fabm%link_horizontal_data(horizontal_id, Ed_0m(jl,:,:)) ! W m^-2

                  varname = 'surf_diffuse_downward_irradiance_' // lam_strings(jl) // '_nm'
                  horizontal_id = model_fabm%get_horizontal_variable_id(TRIM(varname))
                  call model_fabm%link_horizontal_data(horizontal_id, Es_0m(jl,:,:)) ! W m^-2

            end do
      end if
        ! Complete initialization and check whether FABM has all dependencies fulfilled
        ! (i.e., whether all required calls to model%link_*_data have been made)

        write(*,*) 'Finalize initialization and check whether FABM has all dependencies fulfilled ...'

! reduce diagnostis

        do jn = 1, size(model_fabm%interior_diagnostic_variables)
            if (diaWR(jn) .EQ. 1) then
                    model_fabm%interior_diagnostic_variables(jn)%save = .TRUE.
            else
                    model_fabm%interior_diagnostic_variables(jn)%save = .FALSE.
            endif
        enddo
 
        do jn = 1, size(model_fabm%horizontal_diagnostic_variables)
            if (diaWR_2d(jn) .EQ. 1) then
                    model_fabm%horizontal_diagnostic_variables(jn)%save = .TRUE.
            else
                    model_fabm%horizontal_diagnostic_variables(jn)%save = .FALSE.
            endif
        enddo

        call model_fabm%start()
        write(*,*) 'done'

        ! Initialize the tracers
        ! This sets the values of arrays sent to model%link_interior_state_data,
        ! in this case those in interior_state.
        do ji = 1, jpi
           do jj = 1, jpj
              call model_fabm%initialize_interior_state(1, jpk, jj, ji)
           end do
        end do
#ifdef key_trc_fabm
       allocate(tra_DIA(jptra_dia))
       do jn=1,size(model_fabm%interior_diagnostic_variables)
            if (model_fabm%interior_diagnostic_variables(jn)%save) then
                  allocate(tra_DIA(jn)%data(jpk,jpj,jpi))
                  tra_DIA(jn)%data(:,:,:) = huge(tra_DIA(jn)%data(1,1,1))
            endif
       end do
      allocate(tra_DIA_2d(jptra_dia_2d))
      do jn=1,size(model_fabm%horizontal_diagnostic_variables)
            if (model_fabm%horizontal_diagnostic_variables(jn)%save) then
                  allocate(tra_DIA_2d(jn)%data(jpj,jpi))
                  tra_DIA_2d(jn)%data(:,:) = huge(tra_DIA_2d(jn)%data(1,1))
            endif
      end do
             allocate(tra_DIA_IO(jptra_dia))
       do jn=1,size(model_fabm%interior_diagnostic_variables)
            if (model_fabm%interior_diagnostic_variables(jn)%save) then
                  allocate(tra_DIA_IO(jn)%data(jpk,jpj,jpi))
                  tra_DIA_IO(jn)%data(:,:,:) = huge(tra_DIA_IO(jn)%data(1,1,1))
            endif
      end do

       allocate(tra_DIA_IO_HIGH(jptra_dia_high))
       do jn=1, jptra_dia_high
            if (model_fabm%interior_diagnostic_variables(highfreq_table_dia(jn))%save) then
                  allocate(tra_DIA_IO_HIGH(jn)%data(jpk,jpj,jpi))
                  tra_DIA_IO_HIGH(jn)%data(:,:,:) = huge(tra_DIA_IO_HIGH(jn)%data(1,1,1))
            endif
      end do

      allocate(tra_DIA_2d_IO(jptra_dia_2d))
      do jn=1,size(model_fabm%horizontal_diagnostic_variables)
            if (model_fabm%horizontal_diagnostic_variables(jn)%save) then
                  allocate(tra_DIA_2d_IO(jn)%data(jpj,jpi))
                  tra_DIA_2d_IO(jn)%data(:,:) = huge(tra_DIA_2d_IO(jn)%data(1,1))
            endif
      end do
      allocate(tra_DIA_2d_IO_HIGH(jptra_dia2d_high))
      do jn=1,jptra_dia2d_high
            if (model_fabm%horizontal_diagnostic_variables(highfreq_table_dia2d(jn))%save) then
                  allocate(tra_DIA_2d_IO_HIGH(jn)%data(jpj,jpi))
                  tra_DIA_2d_IO_HIGH(jn)%data(:,:) = huge(tra_DIA_2d_IO_HIGH(jn)%data(1,1))
            endif
      end do
#endif
! At this point, initialization is complete.
#endif

#ifdef Mem_Monitor
       mem_all=get_mem(err) - aux_mem
#endif

      END subroutine myalloc_BIO

#ifdef key_trc_fabm

 logical function is_FABM_instance(model,instance_name)
        class (type_fabm_model), pointer, intent(in) :: model
        character(len=*), intent(in) :: instance_name
        type (type_model_list_node), pointer :: instance
   
        instance => model%root%children%first
        
      do while (associated(instance))
         if (instance%model%user_created) then
            if (trim(instance%model%name) == trim(instance_name)) then
                  is_FABM_instance = .true.
                  return
            endif
        end if
        instance => instance%next
      end do

      is_FABM_instance = .false.
      return

      end function is_FABM_instance
#endif

      subroutine clean_memory_bio()

            deallocate(bfm_trn)
            deallocate(bfm_tra)
            deallocate(surf_mask)
            deallocate(ogstm_co2)
            deallocate(co2_IO)
            deallocate(ogstm_sedipi)
            deallocate(ogstm_ph)
            deallocate(NPPF2)

      end subroutine clean_memory_bio



      END MODULE 
