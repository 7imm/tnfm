module tfm_num
  use tfm_essentials 
  use tfm_constants
  use tfm_liquid
  use tfm_temperature
  use tfm_grain
  use tfm_density
  use tfm_llStructure
  implicit none


  ! interface for liquid water
  interface
    subroutine liquid_inter(nz, dt, depth, density, temperature,  &
      & grain_radius, water_content, liquid_accumulation, runoff, &
      & van_genuchten_model)
      use tfm_essentials
      use tfm_liquid
      implicit none

      integer, intent(in)                   :: nz
      real(prec), intent(in)                :: dt
      real(prec), dimension(nz), intent(in) :: depth
      real(prec), dimension(nz), intent(in) :: grain_radius
      real(prec), intent(in)                :: liquid_accumulation 

      procedure(van_genuchten_inter), pointer :: van_genuchten_model

      real(prec), dimension(nz), intent(inout) :: density
      real(prec), dimension(nz), intent(inout) :: temperature
      real(prec), dimension(nz), intent(inout) :: water_content 
      real(prec), intent(inout)                :: runoff
    end subroutine liquid_inter
  end interface


  ! interface for densification function
  interface
    function density_inter(nz, dz, depth, density, temperature, age, grain_radius)
      use tfm_essentials
      implicit none

      integer, intent(in)                   :: nz
      real(prec), intent(in)                :: dz
      real(prec), dimension(nz), intent(in) :: depth
      real(prec), dimension(nz), intent(in) :: density
      real(prec), dimension(nz), intent(in) :: temperature
      real(prec), dimension(nz), intent(in) :: age
      real(prec), dimension(nz), intent(in) :: grain_radius

      real(prec), dimension(nz) :: density_inter
    end function density_inter
  end interface


  ! interface for temperature
  interface
    function temperature_inter(nz, dt, depth, density, temperature, &
      & heat_capacity, thermal_conductivity)
      use tfm_essentials
      implicit none

      integer, intent(in)                   :: nz
      real(prec), intent(in)                :: dt
      real(prec), dimension(nz), intent(in) :: depth
      real(prec), dimension(nz), intent(in) :: density
      real(prec), dimension(nz), intent(in) :: temperature 
      real(prec), dimension(nz), intent(in) :: heat_capacity
      real(prec), dimension(nz), intent(in) :: thermal_conductivity

      real(prec), dimension(nz) :: temperature_inter
    end function temperature_inter
  end interface


  ! interface for heat capacity
  interface
    function heatcap_inter(nz, density, temperature, liquid_water)
      use tfm_essentials
      implicit none
      
      integer, intent(in)                   :: nz
      real(prec), dimension(nz), intent(in) :: density
      real(prec), dimension(nz), intent(in) :: temperature
      real(prec), dimension(nz), intent(in) :: liquid_water
      real(prec), dimension(nz)             :: heatcap_inter
    end function heatcap_inter
  end interface


  ! interface for dry firn thermal conductivity
  interface
    function thermcond_inter(nz, density, temperature)
      use tfm_essentials
      implicit none
      
      integer, intent(in)                   :: nz
      real(prec), dimension(nz), intent(in) :: density
      real(prec), dimension(nz), intent(in) :: temperature
      real(prec), dimension(nz)             :: thermcond_inter
    end function thermcond_inter
  end interface


  ! interface for saturated thermal conductivity
  interface
    function saturation_thermcond_inter(nz, density)
      use tfm_essentials
      implicit none

      integer, intent(in)                   :: nz
      real(prec), dimension(nz), intent(in) :: density
      real(prec), dimension(nz)             :: saturation_thermcond_inter
    end function saturation_thermcond_inter
  end interface


  ! interface for liquid thermal conductivity
  interface
    function liquid_thermcond_inter(nz, density, temperature, &
      & liquid_water, thermcond_model, sat_thermcond_model)
      use tfm_essentials
      implicit none

      integer, intent(in)                            :: nz
      real(prec), dimension(nz), intent(in)          :: density
      real(prec), dimension(nz), intent(in)          :: temperature
      real(prec), dimension(nz), intent(in)          :: liquid_water
      procedure(thermcond_inter), pointer            :: thermcond_model
      procedure(saturation_thermcond_inter), pointer :: sat_thermcond_model

      real(prec), dimension(nz) :: liquid_thermcond_inter
    end function liquid_thermcond_inter
  end interface


  ! interface for grain growth
  interface
    function grain_growth_inter(nz, dt, temperature, density, &
      & liquid_water, grain_radius)
      use tfm_essentials
      implicit none

      integer, intent(in)                   :: nz
      real(prec), intent(in)                :: dt
      real(prec), dimension(nz), intent(in) :: temperature
      real(prec), dimension(nz), intent(in) :: density
      real(prec), dimension(nz), intent(in) :: liquid_water
      real(prec), dimension(nz), intent(in) :: grain_radius

      real(prec), dimension(nz) :: grain_growth_inter
    end function grain_growth_inter
  end interface


  type sim_models
    procedure(density_inter),              pointer, nopass :: dens_model             => null()
    procedure(temperature_inter),          pointer, nopass :: temp_model             => null()
    procedure(heatcap_inter),              pointer, nopass :: heatcap_model          => null()
    procedure(thermcond_inter),            pointer, nopass :: thermcond_model        => null()
    procedure(liquid_thermcond_inter),     pointer, nopass :: liquid_thermcond_model => null()
    procedure(saturation_thermcond_inter), pointer, nopass :: sat_thermcond_model    => null()
    procedure(liquid_inter),               pointer, nopass :: liquid_model           => null()
    procedure(grain_growth_inter),         pointer, nopass :: grain_model            => null()
    procedure(van_genuchten_inter),        pointer, nopass :: van_genuchten_model    => null()
  end type sim_models

  
  contains


  subroutine tfm_num_modelinit(solve_density, solve_temperature, &
    & solve_heat_capacity, solve_thermal_conductivity, &
    & solve_liquid_thermal_conductivity, solve_saturation_thermal_conductivity, &
    & solve_liquid, solve_van_genuchten, solve_grain_growth, models)

    implicit none

    character(len=*), intent(in), optional :: solve_density
    character(len=*), intent(in), optional :: solve_temperature
    character(len=*), intent(in), optional :: solve_heat_capacity
    character(len=*), intent(in), optional :: solve_thermal_conductivity
    character(len=*), intent(in), optional :: solve_liquid_thermal_conductivity
    character(len=*), intent(in), optional :: solve_saturation_thermal_conductivity
    character(len=*), intent(in), optional :: solve_liquid
    character(len=*), intent(in), optional :: solve_van_genuchten
    character(len=*), intent(in), optional :: solve_grain_growth

    ! pointer definition
    type(sim_models), intent(inout) :: models
    
    ! default / fallback
    models%dens_model             => null()
    models%temp_model             => null()
    models%grain_model            => null()
    models%heatcap_model          => null()
    models%thermcond_model        => null()
    models%liquid_thermcond_model => null()
    models%sat_thermcond_model    => null()
    models%liquid_model           => null()
    models%van_genuchten_model    => null()

    if ( present(solve_density) ) then
      if ( solve_density == 'false' ) then
        models%dens_model => null()
      else if ( solve_density == 'medley2020' ) then
        models%dens_model => tfm_density_medley2020
      else if ( solve_density == 'herron1980' ) then
        models%dens_model => tfm_density_herron1980
      else if ( solve_density == 'arthern2010' ) then
        models%dens_model => tfm_density_arthern2010
      else if ( solve_density == 'ligtenberg2011' ) then
        models%dens_model => tfm_density_ligtenberg2011
      else if ( solve_density == 'simonsen2013' ) then
        models%dens_model => tfm_density_simonsen2013
      else if ( solve_density == 'arthern1998' ) then
        models%dens_model => tfm_density_arthern1998
      else if ( solve_density == 'li2003' ) then
        models%dens_model => tfm_density_li2003
      else if ( solve_density == 'helsen2008' ) then
        models%dens_model => tfm_density_helsen2008
      else if ( solve_density == 'breant2017' ) then
        models%dens_model => tfm_density_breant2017
      else if ( solve_density == 'zwinger2007' ) then
        models%dens_model => tfm_density_zwinger2007
      else if ( solve_density == 'greve2009' ) then
        models%dens_model => tfm_density_greve2009
      else if ( solve_density == 'gagliardini1998' ) then
        models%dens_model => tfm_density_gagliardini1998
      else if ( solve_density == 'timmsfit' ) then
        models%dens_model => tfm_density_timmsfit
      else if ( solve_density == 'sintering' ) then
        models%dens_model => tfm_density_sintering
      else
        print *, 'module: tfm_num'
        print *, 'subroutine: tfm_num_modelinit'
        print *, 'can not find densification model: ', solve_density
        print *, 'stopping right here!'
        STOP
      end if
    end if

    if ( present(solve_temperature) ) then
      if ( solve_temperature == 'false' ) then
        models%temp_model => null()
      else if ( solve_temperature == 'true' ) then
        models%temp_model => tfm_temperature_diffusion
      else
        print *, 'module: tfm_num'
        print *, 'subroutine: tfm_num_modelinit'
        print *, 'can not find temperature model: ', solve_temperature
        print *, 'stopping right here!'
        STOP
      end if
    end if
    
    if ( present(solve_liquid) ) then
      if ( solve_liquid == 'false' ) then
        models%liquid_model => null()
      else if ( solve_liquid == 'bucket' ) then
        models%liquid_model => tfm_liquid_bucket
      else if ( solve_liquid == 'richards_equation' ) then
        models%liquid_model => tfm_liquid_RichardsEquation

        if ( present(solve_van_genuchten) ) then
          if ( solve_van_genuchten == 'false' ) then
            print *, 'module: tfm_num'
            print *, 'subroutine: tfm_num_modelinit'
            print *, 'method "richards_equation" is defined for'
            print *, 'solving liquid water, but no model for the'
            print *, 'van Genuchten parameters is given!'
            print *, 'stopping right here!'
            STOP
          else if ( solve_van_genuchten == 'daanen2009' ) then
            models%van_genuchten_model => vgParametersDaanen2009
          else if ( solve_van_genuchten == 'yamaguchi2010' ) then
            models%van_genuchten_model => vgParametersYamaguchi2010
          else if ( solve_van_genuchten == 'yamaguchi2012' ) then
            models%van_genuchten_model => vgParametersYamaguchi2012
          else
            print *, 'module: tfm_num'
            print *, 'subroutine: tfm_num_modelinit'
            print *, 'can not find van Genuchten model:', solve_van_genuchten
            print *, 'stopping right here!'
            STOP
          end if
        
        else
          print *, 'module: tfm_num'
          print *, 'subroutine: tfm_num_modelinit'
          print *, 'methods "richards_equation" is defined for'
          print *, 'solving liquid water, but no model for the'
          print *, 'van Genuchten parameters is given!'
          print *, 'stopping right here!'
          STOP
        end if

      else
        print *, 'module: tfm_num'
        print *, 'subroutine: tfm_num_modelinit'
        print *, 'can not find liquid model: ', solve_liquid
        print *, 'stopping right here!'
        STOP
      end if
    end if

    if ( present(solve_heat_capacity) ) then
      if ( solve_heat_capacity == 'false' ) then
        models%heatcap_model => null()
      else if ( solve_heat_capacity == 'paterson1994' ) then
        models%heatcap_model => tfm_temperature_capacity_paterson1994
      else if ( solve_heat_capacity == 'cuffey2010' ) then
        models%heatcap_model => tfm_temperature_capacity_Cuffey2010
      else
        print *, 'module: tfm_num'
        print *, 'subroutine: tfm_num_modelinit'
        print *, 'can not find heat capacity model: ', solve_heat_capacity
        print *, 'stopping right here!'
        STOP
      end if
    end if

    if ( present(solve_thermal_conductivity) ) then
      if ( solve_thermal_conductivity == 'false' ) then
        models%thermcond_model => null()
      else if ( solve_thermal_conductivity == 'sturm1997' ) then
        models%thermcond_model => tfm_temperature_conduct_sturm1997
      else if ( solve_thermal_conductivity == 'calonne2019' ) then
        models%thermcond_model => tfm_temperature_conduct_calonne2019
      else if ( solve_thermal_conductivity == 'marchenko2019' ) then
        models%thermcond_model => tfm_temperature_conduct_marchenko2019
      else if ( solve_thermal_conductivity == 'miller1969upperbound' ) then
        models%thermcond_model => tfm_temperature_conduct_miller1969UpperBound
      else if ( solve_thermal_conductivity == 'miller1969lowerbound' ) then
        models%thermcond_model => tfm_temperature_conduct_miller1969LowerBound
      else if ( solve_thermal_conductivity == 'geometricmean' ) then
        models%thermcond_model => tfm_temperature_conduct_geomMean
      else
        print *, 'module: tfm_num'
        print *, 'subroutine: tfm_num_modelinit'
        print *, 'can not find thermal conductivity model: ', solve_heat_capacity
        print *, 'stopping right here!'
        STOP
      end if
    end if

    if ( present(solve_liquid_thermal_conductivity) ) then
      if ( solve_liquid_thermal_conductivity == 'false' ) then
        models%liquid_thermcond_model => null()
      else if ( solve_liquid_thermal_conductivity == 'geometricmean' ) then
        models%liquid_thermcond_model => tfm_temperature_liquid_cond_geomMean
      else if ( solve_liquid_thermal_conductivity == 'voigt' ) then
        models%liquid_thermcond_model => tfm_temperature_liquid_cond_voigt
      else
        print *, 'module: tfm_num'
        print *, 'subroutine: tfm_num_modelinit'
        print *, 'can not find liquid thermal conductivity mode: ', solve_liquid_thermal_conductivity
        print *, 'stopping right here!'
        STOP
      end if

      if (                                                    &
      &  present(solve_saturation_thermal_conductivity)       &
      &  .and. (solve_liquid_thermal_conductivity /= 'false') &
      ) then

        if ( solve_saturation_thermal_conductivity == 'geometricmean' ) then
          models%sat_thermcond_model => tfm_temperature_sat_cond_geomMean
        else if ( solve_saturation_thermal_conductivity == 'voigt' ) then
          models%sat_thermcond_model => tfm_temperature_sat_cond_Voigt
        else if ( solve_saturation_thermal_conductivity == 'reuss' ) then
          models%sat_thermcond_model => tfm_temperature_sat_cond_Reuss
        else if ( solve_saturation_thermal_conductivity == 'miller1969upperbound' ) then
          models%sat_thermcond_model => tfm_temperature_sat_cond_Miller1969UpperBound
        else if ( solve_saturation_thermal_conductivity == 'miller1969lowerbound' ) then
          models%sat_thermcond_model => tfm_temperature_sat_cond_Miller1969LowerBound
        else
          print *, 'module: tfm_num'
          print *, 'subroutine: tfm_num_modelinit'
          print *, 'There is a model defined for solving the thermal'
          print *, 'conductivity in presence of liquid water. But '
          print *, 'there is no model defined how to solve the thermal'
          print *, 'conductivity at water saturation.'
          print *, 'Stopping right here!'
          STOP
        end if
      end if
    end if

    if ( present(solve_grain_growth) ) then
      if ( solve_grain_growth == 'false' ) then
        models%grain_model => null()
      else if ( solve_grain_growth == 'arthern2010' ) then
        models%grain_model => tfm_grain_arthern2010
      else if ( solve_grain_growth == 'zwally2002' ) then
        models%grain_model => tfm_grain_zwally2002
      else if ( solve_grain_growth == 'brun1989' ) then
        models%grain_model => tfm_grain_brun1989
      else if ( solve_grain_growth == 'tusima1978' ) then
        models%grain_model => tfm_grain_tusima1978
      else if ( solve_grain_growth == 'katsushima2009' ) then
        models%grain_model => tfm_grain_katsushima2009
      else
        print *, 'module: tfm_num'
        print *, 'subroutine: tfm_num_modelinit'
        print *, 'can not find grain growth model: ', solve_grain_growth
        print *, 'stopping right here!'
        STOP
      end if
    end if
  end subroutine tfm_num_modelinit

  
  subroutine tfm_num_step(dt, models, props, runoff, liquid_acc)
    implicit none

    real(prec), intent(in)              :: dt
    type(sim_models), intent(in)        :: models
    type(llProps), intent(inout)        :: props
    real(prec), intent(inout), optional :: runoff
    real(prec), intent(in), optional    :: liquid_acc

    integer                   :: nz
    integer                   :: n, m
    real(prec), dimension(6)  :: residuum

    real(prec), dimension(props%depth%length)        :: depth
    real(prec), dimension(props%density%length)      :: density
    real(prec), dimension(props%temperature%length)  :: temperature
    real(prec), dimension(props%grain_radius%length) :: grain_radius
    real(prec), dimension(props%heatcap%length)      :: heat_capacity
    real(prec), dimension(props%thermcond%length)    :: thermal_conductivity
    real(prec), dimension(props%liquidwater%length)  :: liquidwater
    real(prec), dimension(props%age%length)          :: age

    real(prec), dimension(props%depth%length)        :: n_depth
    real(prec), dimension(props%density%length)      :: n_density
    real(prec), dimension(props%temperature%length)  :: n_temperature
    real(prec), dimension(props%grain_radius%length) :: n_grain_radius
    real(prec), dimension(props%heatcap%length)      :: n_heat_capacity
    real(prec), dimension(props%thermcond%length)    :: n_thermal_conductivity
    real(prec), dimension(props%liquidwater%length)  :: n_liquidwater
    real(prec), dimension(props%age%length)          :: n_age

    real(prec), dimension(props%depth%length)        :: d_depth
    real(prec), dimension(props%density%length)      :: d_density
    real(prec), dimension(props%temperature%length)  :: d_temperature
    real(prec), dimension(props%grain_radius%length) :: d_grain_radius
    real(prec), dimension(props%density%length)      :: dens_residuum

    ! initialization
    nz                   = props%depth%length
    depth                = llGetData(props%depth)
    density              = llGetData(props%density)
    temperature          = llGetData(props%temperature)
    grain_radius         = llGetData(props%grain_radius)
    heat_capacity        = llGetData(props%heatcap)
    thermal_conductivity = llGetData(props%thermcond)
    liquidwater          = llGetDAta(props%liquidwater)
    age                  = llGetData(props%age)

    ! liquid model
    !if ( (associated(models%liquid_model)) .and. (liquid_acc > 0.0) ) then
    if ( associated(models%liquid_model) ) then
      call models%liquid_model(     &
      &  nz,                        &
      &  dt,                        &
      &  depth,                     &
      &  density,                   &
      &  temperature,               &
      &  grain_radius,              &
      &  liquidwater,               &
      &  liquid_acc,                &
      &  runoff,                    &
      &  models%van_genuchten_model &
      )
    end if

    n_depth                = depth
    n_density              = density
    n_temperature          = temperature
    n_grain_radius         = grain_radius
    n_heat_capacity        = heat_capacity
    n_thermal_conductivity = thermal_conductivity
    n_liquidwater          = liquidwater
    n_age                  = age

    d_density      = 0.0_prec
    d_temperature  = 0.0_prec
    d_grain_radius = 0.0_prec
    residuum       = 0.0_prec
    residuum(1)    = -9999.9_prec

    ! Picard loop
    n = 0
    do while ( (maxval(abs(residuum)) > 1.0e-2_prec) .and. (n < 100) )

      ! density model
      if ( associated(models%dens_model) ) then
        d_density = models%dens_model( &
        &  nz, dt,                     &
        &  depth=n_depth,              &
        &  density=n_density,          &
        &  temperature=n_temperature,  &
        &  age=n_age,                  &
        &  grain_radius=grain_radius   &
        )

        ! There is the possibility that the residuum of the density is
        ! always high because despite the density is converging. This
        ! happens due to the discontinuous function describing
        ! densification. The density "flickers" around the value of
        ! 550 kg m-3. Therefore the residuum at this density is forced
        ! to zero, (which is not ideal).
        dens_residuum = abs(n_density - (density + d_density))
        do m = 1, nz, 1
          if ( floor(n_density(m)) == 550.0_prec       ) dens_residuum(m) = 0.0_prec
          if ( floor(n_density(m)) == CLOSEOFF_DENSITY ) dens_residuum(m) = 0.0_prec
        end do
        residuum(1) = maxval(dens_residuum)
      end if

      ! heat capacity model
      if ( associated(models%heatcap_model) ) then
        n_heat_capacity = models%heatcap_model( &
        &  nz,                                  &
        &  n_density,                           &
        &  n_temperature,                       &
        &  n_liquidwater                        &
        )
        !residuum(2) = maxval(abs(n_heat_capacity - heatcap))
      end if

      ! thermal conductivity model
      if ( associated(models%thermcond_model) ) then
        if (                                                 &
        &  (any(n_liquidwater > 0.0_prec))                   &
        &  .and. (associated(models%liquid_thermcond_model)) &
        ) then
          n_thermal_conductivity = models%liquid_thermcond_model( &
          &  nz,                                                  &
          &  n_density,                                           &
          &  n_temperature,                                       &
          &  n_liquidwater,                                       &
          &  models%thermcond_model,                              &
          &  models%sat_thermcond_model                           &
          )
        else
          n_thermal_conductivity = models%thermcond_model( &
          &  nz,                                           &
          &  n_density,                                    &
          &  n_temperature                                 &
          )
        end if
        !residuum(3) = maxval(abs(n_thermal_conductivity - thermcond))
      end if

      ! temperature model
      if ( associated(models%temp_model) ) then
        d_temperature = models%temp_model(             &
        &  nz, dt,                                     &
        &  depth=n_depth,                              &
        &  density=n_density,                          &
        &  temperature=temperature,                    &
        &  heat_capacity=n_heat_capacity,              &
        &  thermal_conductivity=n_thermal_conductivity &
        )
        residuum(4) = maxval(abs(                        &
        &  n_temperature - (temperature + d_temperature) &
        ))
      end if

      ! grain growth model
      if ( associated(models%grain_model) ) then
        d_grain_radius = models%grain_model( &
        &  nz, dt,                           &
        &  temperature=n_temperature,        &
        &  density=n_density,                &
        &  liquid_water=n_liquidwater,       &
        &  grain_radius=grain_radius         &
        )
        residuum(5) = maxval(abs(                           &
        &  n_grain_radius - (grain_radius + d_grain_radius) &
        ))
      end if

      ! depth evolution
      d_depth = tfm_density_depth( &
      &  nz,                       &
      &  depth=depth,              &
      &  density=density,          &
      &  d_density=d_density       &
      )
      residuum(6) = maxval(abs(n_depth - (depth + d_depth)))

      ! reassignment
      n_depth                = (depth + d_depth)
      n_density              = (density + d_density)
      n_temperature          = (temperature + d_temperature)
      n_grain_radius         = (grain_radius + d_grain_radius)
      n_heat_capacity        = n_heat_capacity
      n_thermal_conductivity = n_thermal_conductivity

      n = n + 1
    end do

    ! raise the age
    n_age = tfm_num_age(nz, dt, age)

    ! new value
    call llUpdateList(props%depth,        n_depth)
    call llUpdateList(props%density,      n_density)
    call llUpdateList(props%temperature,  n_temperature)
    call llUpdateList(props%grain_radius, n_grain_radius)
    call llUpdateList(props%heatcap,      n_heat_capacity)
    call llUpdateList(props%thermcond,    n_thermal_conductivity)
    call llUpdateList(props%age,          n_age)
    call llUpdateList(props%liquidwater,  n_liquidwater)
  end subroutine tfm_num_step


  subroutine tfm_num_surface(dt, forcing, models, props)
    implicit none

    real(prec), intent(in)               :: dt
    real(prec), dimension(7), intent(in) :: forcing
    type(sim_models), intent(in)         :: models

    type(llProps), intent(inout) :: props

    real(prec), dimension(props%depth%length)        :: depth
    real(prec), dimension(props%density%length)      :: density
    real(prec), dimension(props%temperature%length)  :: temperature
    real(prec), dimension(props%heatcap%length)      :: heatcap
    real(prec), dimension(props%thermcond%length)    :: thermcond
    real(prec), dimension(props%grain_radius%length) :: grain_radius
    real(prec), dimension(props%liquidwater%length)  :: liquidwater
    real(prec), dimension(props%age%length)          :: age

    integer    :: nz
    integer    :: n
    real(prec) :: dz, dm, am
    real(prec) :: surf_dens
    real(prec) :: surf_temp
    real(prec) :: surf_grain
    real(prec) :: solid_acc

    nz         = props%depth%length
    surf_temp  = forcing(3)
    surf_dens  = forcing(4)
    solid_acc  = forcing(5)
    surf_grain = forcing(7)

    ! mass to be removed or added
    dm = solid_acc * dt * WATER_DENSITY

    ! the accumulation is zero
    if ( dm == 0.0_prec ) then
      RETURN

    ! theres accumulation and there are still elements available
    else if ( dm > 0.0_prec ) then

      ! height change computed from surface density
      dz = dm / surf_dens

      ! with very small accumulation there can occure precison problem
      ! causing two layers with the same depth
      if ( ((llGetLast(props%depth) + dz) - llGetLast(props%depth)) <= 1.0E-10_prec ) then
      !if ( ((llGetLast(props%depth) + dz) - llGetLast(props%depth)) == 0.0_prec ) then
        RETURN
      end if

      call llAppendData(             &
      &  props%depth,                &
      &  1, (/ llGetLast(props%depth) + dz /) &
      )
      call llAppendData(    &
      &  props%density,     &
      &  1, (/ surf_dens /) &
      )
      call llAppendData(      &
      &  props%temperature,   &
      &  1, (/ surf_temp /)   &
      )
      call llAppendData(       &
      &  props%grain_radius,   &
      &  1, (/ surf_grain /)   &
      )
      call llAppendData(                       &
      &  props%heatcap,                        &
      &  1, (/ models%heatcap_model(           &
      &    1,                                  &
      &    (/ llGetLast(props%density) /),     &
      &    (/ llGetLast(props%temperature) /), &
      &     (/ llGetLast(props%liquidwater) /) &
      &  ) /)                                  &
      )
      call llAppendData(                                                          &
      &  props%thermcond,                                                         &
      &  1, (/ models%thermcond_model(                                            &
      &     1, (/ llGetLast(props%density) /), (/ llGetLast(props%temperature) /) &
      &   ) /)                                                                    &
      )
      call llAppendData(      &
      &  props%liquidwater,   &
      &  1, (/ 0.0_prec /)    &
      )
      call llAppendData(   &
      &  props%age,        &
      &  1, (/ 0.0_prec /) &
      )

    ! theres ablation
    else if ( dm < 0.0_prec ) then
      
      depth        = llGetData(props%depth)
      density      = llGetData(props%density)
      temperature  = llGetData(props%temperature)
      heatcap      = llGetData(props%heatcap)
      thermcond    = llGetData(props%thermcond)
      grain_radius = llGetData(props%grain_radius)
      liquidwater  = llGetData(props%liquidwater)
      age          = llGetData(props%age)

      ! removal of layers
      am = -dm
      do n = nz - 1, 1, -1
        am = am - ((depth(n+1) - depth(n)) * density(n))
        if ( am <= 0.0_prec ) EXIT
        dm = dm + ((depth(n+1) - depth(n)) * density(n))
      end do
      dz = dm / density(n)
      n = n + 1

      ! interpolation
      call tfm_num_lin_interp(               &
      &  depth(n), depth(n-1),               &
      &  temperature(n), temperature(n-1),   &
      &  dz, temperature(n)                  &
      )
      call tfm_num_lin_interp(               &
      &  depth(n), depth(n-1),               &
      &  grain_radius(n), grain_radius(n-1), &
      &  dz, grain_radius(n)                 &
      )
      call tfm_num_lin_interp(               &
      &  depth(n), depth(n-1),               &
      &  liquidwater(n), liquidwater(n-1),   &
      &  dz, liquidwater(n)                  &
      )
      call tfm_num_lin_interp(               &
      &  depth(n), depth(n-1),               &
      &  heatcap(n), heatcap(n-1),           &
      &  dz, heatcap(n)                      &
      )
      call tfm_num_lin_interp(               &
      &  depth(n), depth(n-1),               &
      &  thermcond(n), thermcond(n-1),       &
      &  dz, thermcond(n)                    &
      )
      call tfm_num_lin_interp(               &
      &  depth(n), depth(n-1),               &
      &  age(n), age(n-1),                   &
      &  dz, age(n)                          &
      )

      ! new value
      call llUpdateList(props%depth,        depth)
      call llUpdateList(props%density,      density)
      call llUpdateList(props%temperature,  temperature)
      call llUpdateList(props%grain_radius, grain_radius)
      call llUpdateList(props%heatcap,      heatcap)
      call llUpdateList(props%thermcond,    thermcond)
      call llUpdateList(props%age,          age)
      call llUpdateList(props%liquidwater,  liquidwater)
 
      ! fill removed layers with NaN value
      call llDropData(props%depth,        -(nz - n))
      call llDropData(props%density,      -(nz - n))
      call llDropData(props%temperature,  -(nz - n))
      call llDropData(props%heatcap,      -(nz - n))
      call llDropData(props%thermcond,    -(nz - n))
      call llDropData(props%grain_radius, -(nz - n))
      call llDropData(props%liquidwater,  -(nz - n))
      call llDropData(props%age,          -(nz - n))

      ! new depth / height of the uppermost layer
      !props%depth%tail%data(props%depth%tind - 1) = (     &
      !&  props%depth%tail%data(props%depth%tind - 1) + dz &
      !)
      depth(n) = (depth(n) + dz)
      call llUpdateList(props%depth, depth)
    end if
  end subroutine tfm_num_surface


  subroutine tfm_num_lin_interp(z0, z1, v0, v1, dz, v)
    implicit none
    
    real(prec), intent(in) :: z0, z1
    real(prec), intent(in) :: v0, v1
    real(prec), intent(in) :: dz

    real(prec), intent(inout) :: v

    v = v0 + ((v1 - v0) / (z1 - z0)) * dz
  end subroutine tfm_num_lin_interp


  function tfm_num_age(nz, dt, age) result (n_age)
    implicit none

    integer, intent(in)                   :: nz
    real(prec), intent(in)                :: dt
    real(prec), dimension(nz), intent(in) :: age
    real(prec), dimension(nz)             :: n_age

    n_age = age + dt
  end function tfm_num_age


  subroutine tfm_num_trimProfileLength(props, length)
    use tfm_llStructure
    implicit none

    type(llProps), intent(inout) :: props
    real(prec), intent(in)       :: length

    integer                                   :: n
    real(prec), dimension(props%depth%length) :: depth

    depth = llGetData(props%depth)
    do n = 1, size(depth), 2
      if ( (depth(size(depth)) - depth(n)) <= length ) EXIT
    end do
    n = (n - 1)

    call llPropsDropData(props, n)
  end subroutine tfm_num_trimProfileLength


  subroutine tfm_num_trimProfileAge(props, max_age)
    use tfm_llStructure
    implicit none

    type(llProps), intent(inout) :: props
    real(prec), intent(in)       :: max_age

    integer                                 :: n
    real(prec), dimension(props%age%length) :: age

    age = (llGetData(props%age) / SECONDS_YEAR)
    do n = 1, size(age), 2
      if ( age(n) <= max_age ) EXIT
    end do
    n = (n - 1)

    call llPropsDropData(props, n)
  end subroutine tfm_num_trimProfileAge
end module tfm_num
