module tfm_num
  use settings
  use tfm_constants
  use tfm_liquid
  use tfm_temperature
  use tfm_density
  implicit none


  ! interface for bucket scheme
  interface
    subroutine liquid_inter(nz, dt, depth, density, temperature, &
      &liquid_water, infiltration_rate, runoff)
      use settings
      implicit none

      integer, intent(in)                   :: nz
      real(prec), intent(in)                :: dt
      real(prec), dimension(nz), intent(in) :: depth
      real(prec), intent(in)                :: infiltration_rate

      real(prec), dimension(nz), intent(inout) :: density
      real(prec), dimension(nz), intent(inout) :: temperature
      real(prec), dimension(nz), intent(inout) :: liquid_water
      real(prec), intent(inout)                :: runoff
    end subroutine liquid_inter
  end interface


  ! interface for densification function
  interface
    function density_inter(nz, dz, accumulation, depth, density, temperature)
      use settings
      implicit none

      integer, intent(in)                   :: nz
      real(prec), intent(in)                :: dz
      real(prec), intent(in)                :: accumulation
      real(prec), dimension(nz), intent(in) :: depth
      real(prec), dimension(nz), intent(in) :: density
      real(prec), dimension(nz), intent(in) :: temperature

      real(prec), dimension(nz) :: density_inter
    end function density_inter
  end interface


  ! interface for temperature
  interface
    function temperature_inter(nz, dt, depth, density, temperature, &
      & heat_capacity, thermal_conductivity)
      use settings
      implicit none

      integer, intent(in) :: nz
      real(prec), intent(in) :: dt
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
    function heatcap_inter(nz)
      use settings
      implicit none
      
      integer, intent(in)       :: nz
      real(prec), dimension(nz) :: heatcap_inter
    end function heatcap_inter
  end interface


  ! interface for thermal conductivity
  interface
    function thermcond_inter(nz, density)
      use settings
      implicit none
      
      integer, intent(in)                   :: nz
      real(prec), dimension(nz), intent(in) :: density
      real(prec), dimension(nz)             :: thermcond_inter
    end function thermcond_inter
  end interface


  type sim_models
    procedure(density_inter),     pointer, nopass :: dens_model      => null()
    procedure(temperature_inter), pointer, nopass :: temp_model      => null()
    procedure(heatcap_inter),     pointer, nopass :: heatcap_model   => null()
    procedure(thermcond_inter),   pointer, nopass :: thermcond_model => null()
    procedure(liquid_inter),      pointer, nopass :: liquid_model    => null()
  end type sim_models

  
  contains


  subroutine tfm_num_modelinit(solve_density, solve_temperature, &
    & solve_heat_capacity, solve_thermal_conductivity, solve_liquid, models)

    implicit none

    character(len=*), intent(in), optional :: solve_density
    character(len=*), intent(in), optional :: solve_temperature
    character(len=*), intent(in), optional :: solve_heat_capacity
    character(len=*), intent(in), optional :: solve_thermal_conductivity
    character(len=*), intent(in), optional :: solve_liquid

    ! pointer definition
    type(sim_models), intent(inout) :: models
    
    ! default / fallback
    models%dens_model      => null()
    models%temp_model      => null()
    models%heatcap_model   => null()
    models%thermcond_model => null()
    models%liquid_model    => null()

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
      else
        print *, 'module: tfm_num'
        print *, 'subroutine: tfm_num_modelinit'
        print *, 'can not find heat capacity model: ', solve_heat_capacity
        print *, 'stopping right here!'
      end if
    end if

    if ( present(solve_thermal_conductivity) ) then
      if ( solve_thermal_conductivity == 'false' ) then
        models%thermcond_model => null()
      else if ( solve_thermal_conductivity == 'sturm2007' ) then
        models%thermcond_model => tfm_temperature_conduct_sturm2007
      else
        print *, 'module: tfm_num'
        print *, 'subroutine: tfm_num_modelinit'
        print *, 'can not find thermal conductivity model: ', solve_heat_capacity
        print *, 'stopping right here!'
      end if
    end if
  end subroutine tfm_num_modelinit

  
  subroutine tfm_num_step(nz, dt, models, props, runoff, liquid_acc, solid_acc)

    implicit none

    integer, intent(in)    :: nz
    real(prec), intent(in) :: dt

    type(sim_models), intent(in)               :: models
    real(prec), dimension(6,nz), intent(inout) :: props
    real(prec), intent(inout), optional        :: runoff

    real(prec), intent(in), optional :: liquid_acc
    real(prec), intent(in), optional :: solid_acc

    type(sim_props)           :: p
    integer                   :: n, m
    real(prec), dimension(4)  :: residuum
    real(prec), dimension(nz) :: dens_residuum
    real(prec), dimension(nz) :: d_density, n_density
    real(prec), dimension(nz) :: d_temperature, n_temperature
    real(prec), dimension(nz) :: n_heat_capacity
    real(prec), dimension(nz) :: n_thermal_conductivity

    call tfm_num_assign(nz, props, p)

    ! initialization
    n_density = p%density
    n_temperature = p%temperature
    n_heat_capacity = p%heatcap
    n_thermal_conductivity = p%thermcond
    d_density = 0.0
    d_temperature = 0.0
    residuum = 0.0
    residuum(1) = -9999.9

    ! liquid model
    if ( associated(models%liquid_model) ) then
      call models%liquid_model(nz, dt, p%depth, p%density, &
        & p%temperature, p%liquidwater, liquid_acc, runoff)
    end if

    ! Picard loop
    n = 0
    do while ( (maxval(abs(residuum)) > 1.0e-1) .and. (n < 10000) )
      
      ! density model
      if ( associated(models%dens_model) ) then
        d_density = models%dens_model(            &
        &  nz, dt,                                &
        &  accumulation=(solid_acc + liquid_acc), &
        &  depth=p%depth,                         &
        &  density=n_density,                     &
        &  temperature=n_temperature              &
        )

        ! There is the possibility that the residuum of the density is
        ! always high because despite the density is converging. This
        ! happens due to the discontinuous function describing
        ! densification. The density "flickers" around the value of
        ! 550 kg m-3. Therefore the residuum at this density is forced
        ! to zero, (which is not ideal).
        dens_residuum = abs(n_density - (p%density + d_density))
        do m = 1, nz, 1
          if ( floor(p%density(m)) == 549 ) dens_residuum(m) = 0.0
        end do
        residuum(1) = maxval(dens_residuum)
      end if

      ! heat capacity model
      if ( associated(models%heatcap_model) ) then
        n_heat_capacity = models%heatcap_model(nz)
        !residuum(2) = maxval(abs(n_heat_capacity - p%heatcap))
      end if

      ! thermal conductivity model
      if ( associated(models%thermcond_model) ) then
        n_thermal_conductivity = models%thermcond_model(nz, n_density)
        !residuum(3) = maxval(abs(n_thermal_conductivity - p%thermcond))
      end if

      ! temperature model
      if ( associated(models%temp_model) ) then
        d_temperature = models%temp_model(             &
        &  nz, dt,                                     &
        &  depth=p%depth,                              &
        &  density=n_density,                          &
        &  temperature=p%temperature,                  &
        &  heat_capacity=n_heat_capacity,              &
        &  thermal_conductivity=n_thermal_conductivity &
        )
        residuum(4) = maxval(abs(                          &
        &  n_temperature - (p%temperature + d_temperature) &
        ))
      end if

      ! reassignment
      n_density              = (p%density + d_density)
      n_temperature          = (p%temperature + d_temperature)
      n_heat_capacity        = n_heat_capacity
      n_thermal_conductivity = n_thermal_conductivity

      n = n + 1
    end do

    ! new value
    p%density     = n_density
    p%temperature = n_temperature
    p%heatcap     = n_heat_capacity
    p%thermcond   = n_thermal_conductivity
  end subroutine tfm_num_step


  subroutine tfm_num_surface(nz, np, dt, forcing, models, props)
    implicit none

    integer, intent(in)                  :: np
    real(prec), intent(in)               :: dt
    real(prec), dimension(6), intent(in) :: forcing
    type(sim_models), intent(in)         :: models

    integer, intent(inout)                     :: nz
    real(prec), dimension(6,np), intent(inout) :: props

    type(sim_props) :: p
    integer         :: n
    real(prec)      :: dz, dm, am
    real(prec)      :: surf_temp, surf_dens, solid_acc


    call tfm_num_assign(np, props, p)

    surf_temp = forcing(3)
    surf_dens = forcing(4)
    solid_acc = forcing(5)

    ! mass to be removed
    dm = solid_acc * dt * WATER_DENSITY

    ! theres accumulation and the maximum number of elements is reached
    if ( nz == np .and. dm > 0.0 ) then

      p%depth(1:nz-1)       = p%depth(2:nz)
      p%density(1:nz-1)     = p%density(2:nz)
      p%temperature(1:nz-1) = p%temperature(2:nz)
      p%liquidwater(1:nz-1) = p%liquidwater(2:nz)

      p%depth(nz)       = -9999.9
      p%density(nz)     = -9999.9
      p%temperature(nz) = -9999.9
      p%liquidwater(nz) = -9999.9
      p%heatcap(nz)     = -9999.9
      p%thermcond       = -9999.9

      nz = nz - 1
    end if

    ! theres accumulation and there are still elements available
    if ( dm > 0.0 ) then

      ! height change computed from surface density
      dz = dm / surf_dens

      ! add new layer
      nz = nz + 1
      p%depth(nz)       = p%depth(nz-1) + dz
      p%density(nz)     = surf_dens
      p%temperature(nz) = surf_temp
      p%heatcap(1:nz)   = models%heatcap_model(nz)
      p%thermcond(1:nz) = models%thermcond_model(nz, p%density(1:nz))
      p%liquidwater(nz) = 0.0

    ! theres ablation
    else if ( dm < 0.0 ) then

      ! removal of layers
      am = -dm
      do n = nz - 1, 1, -1
        am = am - ((p%depth(n+1) - p%depth(n)) * p%density(n))
        if ( am <= 0.0 ) EXIT
        dm = dm + ((p%depth(n+1) - p%depth(n)) * p%density(n))
      end do
      dz = dm / p%density(n)
      n = n + 1

      ! interpolation
      call tfm_num_lin_interp(                 &
      &  p%depth(n), p%depth(n-1),             &
      &  p%temperature(n), p%temperature(n-1), &
      &  dz, p%temperature(n)                  &
      )
      call tfm_num_lin_interp(                 &
      &  p%depth(n), p%depth(n-1),             &
      &  p%liquidwater(n), p%liquidwater(n-1), &
      &  dz, p%liquidwater(n)                  &
      )
      call tfm_num_lin_interp(                 &
      &  p%depth(n), p%depth(n-1),             &
      &  p%heatcap(n), p%heatcap(n-1),         &
      &  dz, p%heatcap(n)                      &
      )
      call tfm_num_lin_interp(                 &
      &  p%depth(n), p%depth(n-1),             &
      &  p%thermcond(n), p%thermcond(n-1),     &
      &  dz, p%thermcond(n)                    &
      )

      ! new depth / height of the uppermost layer
      p%depth(n) = p%depth(n) + dz
      nz = n

      ! fill removed layers with NaN value
      p%depth(nz+1:np)       = -9999.9
      p%density(nz+1:np)     = -9999.9
      p%temperature(nz+1:np) = -9999.9
      p%liquidwater(nz+1:np) = -9999.9
      p%heatcap(nz+1:np)     = -9999.9
      p%thermcond(nz+1:np)   = -9999.9
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


  subroutine tfm_num_assign(nz, props, p)
    implicit none

    integer, intent(in)                             :: nz
    real(prec), dimension(6,nz), intent(in), target :: props
    type(sim_props), intent(inout)                  :: p

    p%depth       => props(1,:)
    p%density     => props(2,:)
    p%temperature => props(3,:)
    p%heatcap     => props(4,:)
    p%thermcond   => props(5,:)
    p%liquidwater => props(6,:)
  end subroutine tfm_num_assign
end module tfm_num
