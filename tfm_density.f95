module tfm_density
  use settings
  use tfm_constants
  implicit none

  contains

  subroutine tfm_density_HLtype(nz, dt, comp_accumulation, stage1_params, &
    & stage2_params, temperature, density, d_density)
    implicit none

    integer, intent(in)                   :: nz
    real(prec), intent(in)                :: dt
    real(prec), intent(in)                :: comp_accumulation
    real(prec), dimension(3), intent(in)  :: stage1_params
    real(prec), dimension(3), intent(in)  :: stage2_params
    real(prec), dimension(nz), intent(in) :: temperature
    real(prec), dimension(nz), intent(in) :: density

    real(prec), dimension(nz), intent(inout) :: d_density

    integer                   :: n
    real(prec), dimension(nz) :: c

    ! boundary between first and seconds stage
    do n = nz, 1, -1
      if ( density(n) >= 550.0 ) EXIT
    end do

    ! first stage
    c(n:nz) = (                                                   &
    &  stage1_params(1)                                           &
    &  * (comp_accumulation**stage1_params(2))                    &
    &  * ACC_GRAVITY                                              &
    &  * exp(-stage1_params(3) / (GAS_CONST * temperature(n:nz))) &
    )

    ! seconds stage
    c(1:n-1) = (                                                   &
    &  stage2_params(1)                                            &
    &  * (comp_accumulation**stage2_params(2))                     &
    &  * ACC_GRAVITY                                               &
    &  * exp(-stage2_params(3) / (GAS_CONST * temperature(1:n-1))) &
    )
    
    ! density change
    d_density = ((dt / SECONDS_YEAR) * (c * (ICE_DENSITY - density)))
  end subroutine tfm_density_HLtype


  function tfm_density_medley2020(nz, dt, accumulation, &
    & depth, density, temperature) result(d_density)
    implicit none

    ! first stage parameters
    real(prec), parameter :: ALPHA0 = 0.9250
    real(prec), parameter :: EC0    = 60000.0
    real(prec)            :: a0

    ! seconds stage parameters
    real(prec), parameter :: ALPHA1 = 0.6354
    real(prec), parameter :: EC1    = 56973.0
    real(prec)            :: a1

    ! other parameters
    real(prec), parameter :: EG = 42400.0

    integer, intent(in)                   :: nz
    real(prec), intent(in)                :: dt
    real(prec), intent(in)                :: accumulation
    real(prec), dimension(nz), intent(in) :: depth
    real(prec), dimension(nz), intent(in) :: density
    real(prec), dimension(nz), intent(in) :: temperature

    real(prec), dimension(nz) :: d_density

    integer    :: n
    real(prec) :: comp_accumulation
    real(prec) :: mean_temperature

    ! 10 m temperature
    do n = nz, 1, -1
      if ( (depth(n) - depth(nz)) <= -10.0 ) EXIT
    end do
    mean_temperature = temperature(n)

    ! factor depending on the mean annual temperature
    a0 = 0.07 * exp(EG / (GAS_CONST * mean_temperature))
    a1 = 0.03 * exp(EG / (GAS_CONST * mean_temperature))

    ! change accumulation to (kg a-1 m-2)
    comp_accumulation = accumulation * SECONDS_YEAR * WATER_DENSITY
    if ( comp_accumulation < 0.0 ) comp_accumulation = 0.0

    ! call Herron & Langway model
    call tfm_density_HLtype(  &
    &  nz, dt,                &
    &  comp_accumulation,     &
    &  (/ a0, ALPHA0, EC0 /), &
    &  (/ a1, ALPHA1, EC1 /), &
    &  temperature,         &
    &  density,             &
    &  d_density              &
    )
  end function tfm_density_medley2020


  function tfm_density_herron1980(nz, dt, accumulation, depth, density, &
    & temperature) result(d_density)
    implicit none

    ! first stage parameters
    real(prec), parameter :: ALPHA0 = 1.0
    real(prec), parameter :: EC0    = 10160.0
    real(prec), parameter :: A0     = 11.0

    ! second stage parameters
    real(prec), parameter :: ALPHA1 = 0.5
    real(prec), parameter :: EC1    = 21400.0
    real(prec), parameter :: A1     = 575.0

    integer, intent(in)                   :: nz
    real(prec), intent(in)                :: dt
    real(prec), intent(in)                :: accumulation
    real(prec), dimension(nz), intent(in) :: depth
    real(prec), dimension(nz), intent(in) :: density
    real(prec), dimension(nz), intent(in) :: temperature

    real(prec), dimension(nz) :: d_density

    real(prec) :: comp_accumulation, dummy

    dummy = depth(1)

    ! change accumulation to (m weq. a-1)
    comp_accumulation = accumulation * SECONDS_YEAR
    if ( comp_accumulation < 0.0 ) comp_accumulation = 0.0

    ! call Herron & Langway model
    call tfm_density_HLtype(  &
    &  nz, dt,                &
    &  comp_accumulation,     &
    &  (/ A0, ALPHA0, EC0 /), &
    &  (/ A1, ALPHA1, EC1 /), &
    &  temperature,           &
    &  density,               &
    &  d_density              &
    )
  end function tfm_density_herron1980


  function tfm_density_arthern2010(nz, dt, accumulation, depth, &
    density, temperature) result(d_density)
    implicit none

    ! first stage parameters
    real(prec), parameter :: ALPHA0 = 1.0
    real(prec), parameter :: EC0    = 60000.0
    real(prec)            :: a0

    ! second stage parameters
    real(prec), parameter :: ALPHA1 = 1.0
    real(prec), parameter :: EC1    = 60000.0
    real(prec)            :: a1

    ! further parameters
    real(prec), parameter :: EG = 42400.0

    integer, intent(in)    :: nz
    real(prec), intent(in) :: dt
    real(prec), intent(in) :: accumulation
    real(prec), dimension(nz), intent(in) :: depth
    real(prec), dimension(nz), intent(in) :: density
    real(prec), dimension(nz), intent(in) :: temperature

    real(prec), dimension(nz) :: d_density

    integer    :: n
    real(prec) :: comp_accumulation
    real(prec) :: mean_temperature

    ! change accumulation to (kg a-1 m-2)
    comp_accumulation = accumulation * SECONDS_YEAR * WATER_DENSITY
    if ( comp_accumulation < 0.0 ) comp_accumulation = 0.0

    ! 10 m temperature
    do n = nz, 1, -1
      if ( (depth(n) - depth(nz)) <= -10.0 ) EXIT
    end do
    mean_temperature = temperature(n)

    ! factor depending on the mean annual temperature
    a0 = 0.07 * exp(EG / (GAS_CONST * mean_temperature))
    a1 = 0.03 * exp(EG / (GAS_CONST * mean_temperature))

    ! call Herron & Langway model
    call tfm_density_HLtype(  &
    &  nz, dt,                &
    &  comp_accumulation,     &
    &  (/ a0, ALPHA0, EC0 /), &
    &  (/ a1, ALPHA1, EC1 /), &
    &  temperature,           &
    &  density,               &
    &  d_density              &
    )
  end function tfm_density_arthern2010


  function tfm_density_ligtenberg2011(nz, dt, accumulation, depth, &
    & density, temperature) result(d_density)
    implicit none

    integer, intent(in)                   :: nz
    real(prec), intent(in)                :: dt
    real(prec), intent(in)                :: accumulation
    real(prec), dimension(nz), intent(in) :: depth
    real(prec), dimension(nz), intent(in) :: density
    real(prec), dimension(nz), intent(in) :: temperature

    real(prec), dimension(nz) :: d_density

    integer    :: n
    real(prec) :: comp_accumulation

    d_density = tfm_density_arthern2010(nz, dt, accumulation, depth, &
      & density, temperature)

    ! boundary between first and seconds stage
    do n = nz, 1, -1
      if ( density(n) >= 550.0 ) EXIT
    end do

    ! change accumulation to (kg a-1 m-2)
    comp_accumulation = accumulation * SECONDS_YEAR * WATER_DENSITY
    if ( comp_accumulation < 0.0 ) comp_accumulation = 0.0

    d_density(n:nz) = (                             &
    &  d_density(n:nz)                              &
    &  * (1.435 - (0.151 * log(comp_accumulation))) &
    )
    d_density(1:n-1) = (                            &
    &  d_density(1:nz)                              &
    &  * (2.366 - (0.293 * log(comp_accumulation))) &
    )
  end function tfm_density_ligtenberg2011


  function tfm_density_simonsen2013(nz, dt, accumulation, depth, &
    & density, temperature) result(d_density)
    implicit none

    ! fist stage parameters (as implemented in CFM)
    real(prec), parameter :: F0 = 0.8

    ! second stage parameters (as implemented in CFm)
    real(prec), parameter :: F1 = 1.25

    integer, intent(in)                   :: nz
    real(prec), intent(in)                :: dt
    real(prec), intent(in)                :: accumulation
    real(prec), dimension(nz), intent(in) :: depth
    real(prec), dimension(nz), intent(in) :: density
    real(prec), dimension(nz), intent(in) :: temperature

    real(prec), dimension(nz) :: d_density

    integer    :: n
    real(prec) :: comp_accumulation
    real(prec) :: mean_temperature

    d_density = tfm_density_arthern2010(nz, dt, accumulation, depth, &
      & density, temperature)

    ! boundary between first and seconds stage
    do n = nz, 1, -1
      if ( density(n) >= 550.0 ) EXIT
    end do

    ! change accumulation to (kg a-1 m-2)
    comp_accumulation = accumulation * SECONDS_YEAR * WATER_DENSITY
    if ( comp_accumulation < 0.0 ) comp_accumulation = 0.0

    ! 10 m temperature
    do n = nz, 1, -1
      if ( (depth(n) - depth(nz)) <= -10.0 ) EXIT
    end do
    mean_temperature = temperature(n)

    d_density(n:nz) = F0 * d_density(n:nz)
    d_density(1:n-1) = (                               &
    &  d_density(1:n-1)                                &
    &  * (61.7 / (comp_accumulation**0.5))             &
    &  * exp(-3800.0 / (GAS_CONST * mean_temperature)) &
    )
  end function tfm_density_simonsen2013


  subroutine tfm_density_herron1980_analytical(nz, temperature, &
    & accumulation, surface_density, depth, density)
    implicit none

    integer, intent(in)                   :: nz
    real(prec), intent(in)                :: temperature
    real(prec), intent(in)                :: accumulation
    real(prec), intent(in)                :: surface_density
    real(prec), dimension(nz), intent(in) :: depth

    real(prec), dimension(nz), intent(inout) :: density

    integer    :: n
    real(prec) :: k0, k1, z550
    real(prec), dimension(nz) :: comp_depth

    comp_depth = -(depth - depth(nz))

    k0 =  11.0 * exp(-10160.0 / (GAS_CONST * temperature))
    k1 = 575.0 * exp(-21400.0 / (GAS_CONST * temperature))

    density = (                                            &
    &  (surface_density / (ICE_DENSITY - surface_density)) &
    &  * exp((ICE_DENSITY / 1000.0) * k0 * comp_depth)     &
    )
    density = ICE_DENSITY * (density / (1.0 + density))

    do n = nz, 1, -1
      if ( density(n) >= 550.0 ) EXIT
    end do

    z550 = (                                                               &
    &  comp_depth(n+1)                                                     &
    &  + ((comp_depth(n) - comp_depth(n+1)) / (density(n) - density(n+1))) &
    &  * (550.0 - density(n+1))                                            &
    )

    density(1:n) = (                            &
    &  (550.0 / (ICE_DENSITY - 550.0))          &
    &  * exp(                                   &
    &    (ICE_DENSITY / 1000.0)                 &
    &    * (k1 / (accumulation * SECONDS_YEAR)) &
    &    * (comp_depth(1:n) - z550)             &
    &  )                                        &
    )
    density(1:n) = ICE_DENSITY * (density(1:n) / (1.0 + density(1:n)))
  end subroutine tfm_density_herron1980_analytical
end module tfm_density
