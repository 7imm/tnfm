module tfm_density_tools
  use settings
  use tfm_constants
  implicit none

  interface tfm_density_arrhenius
    module procedure tfm_density_arrhenius_i, tfm_density_arrhenius_n
  end interface tfm_density_arrhenius

  contains


  function tfm_density_arrhenius_i(i, factor, activation_energy, temperature) &
    & result(arrhenius)
    implicit none
    
    integer, intent(in)    :: i
    real(prec), intent(in) :: factor
    real(prec), intent(in) :: activation_energy
    real(prec), intent(in) :: temperature

    integer    :: m
    real(prec) :: arrhenius

    m = 1 * i

    arrhenius = factor * exp(-activation_energy / (GAS_CONST * temperature))
  end function tfm_density_arrhenius_i


  function tfm_density_arrhenius_n(n, factor, activation_energy, temperature) &
    & result(arrhenius)
    implicit none

    integer, intent(in)                  :: n
    real(prec), intent(in)               :: factor
    real(prec), intent(in)               :: activation_energy
    real(prec), dimension(n), intent(in) :: temperature

    real(prec), dimension(n) :: arrhenius

    arrhenius = factor * exp(-activation_energy / (GAS_CONST * temperature))
  end function tfm_density_arrhenius_n


  subroutine tfm_density_lin_interp(z0, z1, v0, v1, dz, v)
    implicit none

    real(prec), intent(in) :: z0, z1
    real(prec), intent(in) :: v0, v1
    real(prec), intent(in) :: dz

    real(prec), intent(inout) :: v

    v = v0 + ((v1 - v0) / (z1 - z0)) * dz
  end subroutine tfm_density_lin_interp


  subroutine tfm_density_do_nothing(nz, variable)
    implicit none

    integer, intent(in) :: nz
    real(prec), dimension(nz), intent(in) :: variable
    real(prec), dimension(nz)             :: nothing

    nothing = variable
  end subroutine tfm_density_do_nothing


  subroutine tfm_density_bagmean(nz, depth, density, dz, mz, bagmean)
    implicit none

    integer, intent(in)                        :: nz
    real(prec), dimension(nz), intent(in)      :: depth
    real(prec), dimension(nz), intent(in)      :: density
    real(prec), intent(in)                     :: dz
    integer, intent(in)                        :: mz
    real(prec), dimension(2,mz), intent(inout) :: bagmean

    integer    :: n
    integer    :: m
    real(prec) :: d

    bagmean = 0.0

    n = nz
    do m = 1, mz, 1
      
      d = 0.0
      bagmean(1,mz-m+1) = depth(nz) - (m * dz) + (0.5 * dz)

      do while ( depth(n) > (depth(nz) - (m * dz)) )
        bagmean(2,mz-m+1) = bagmean(2,mz-m+1) + density(n)
        n = n - 1
        d = d + 1.0
      end do

      bagmean(2,mz-m+1) = bagmean(2,mz-m+1) / d
    end do
  end subroutine tfm_density_bagmean
end module tfm_density_tools



module tfm_density_stress
  use settings
  use tfm_constants
  implicit none

  contains


  subroutine tfm_density_computeStress(nz, depth, density, stress)
    implicit none

    integer, intent(in)                      :: nz
    real(prec), dimension(nz), intent(in)    :: depth
    real(prec), dimension(nz), intent(in)    :: density
    real(prec), dimension(nz), intent(inout) :: stress

    integer                   :: n
    real(prec), dimension(nz) :: dz

    dz(nz) = 0.0_prec
    dz(1:nz-1) =  depth(2:nz) - depth(1:nz-1)

    stress = (dz * density * ACC_GRAVITY)
    do n = nz - 1, 1, -1
      stress(n) = stress(n) + stress(n+1)
    end do
  end subroutine tfm_density_computeStress


  function tfm_density_boyleMariotte(density) result(boyle_mariotte)
    implicit none
    
    real(prec), intent(in) :: density
    real(prec)             :: boyle_mariotte

    boyle_mariotte = (                                &
    &  (density * (ICE_DENSITY - CLOSEOFF_DENSITY))   &
    &  / (CLOSEOFF_DENSITY * (ICE_DENSITY - density)) &
    )
  end function tfm_density_boyleMariotte


  function tfm_density_rel_boyleMariotte(rel_density) result(rel_boyle_mariotte)
    implicit none

    real(prec), intent(in) :: rel_density
    real(prec)             :: rel_boyle_mariotte

    rel_boyle_mariotte = (                                             &
    &  rel_density * (1.0_prec - (CLOSEOFF_DENSITY / ICE_DENSITY))     &
    &  / ((CLOSEOFF_DENSITY / ICE_DENSITY) * (1.0_prec - rel_density)) &
    )
  end function tfm_density_rel_boyleMariotte
end module tfm_density_stress




module tfm_density_processes
  use settings
  use tfm_constants
  use tfm_density_tools
  use tfm_density_stress
  implicit none

  ! parameters
  real(prec), parameter :: DELTA_B      = 9.0e-10_prec
  real(prec), parameter :: OMEGA        = 3.27e-29_prec
  real(prec), parameter :: AIR_PRESSURE = 101325.0_prec
  real(prec), parameter :: H            = 4.0e-6_prec
  real(prec), parameter :: MU           = 0.7_prec

  contains


  function tfm_density_grainBoundarySliding(density, temperature, &
    & grain_radius, stress) result(strain_rate)
    implicit none

    real(prec), intent(in) :: density
    real(prec), intent(in) :: temperature
    real(prec), intent(in) :: grain_radius
    real(prec), intent(in) :: stress

    real(prec) :: strain_rate
    real(prec) :: arrhenius

    ! paramters
    real(prec) :: ABDC    = 3.0e-2_prec
    real(prec) :: QBDC    = 44100.0_prec

    strain_rate = 0.0_prec

    if ( ( density > 0.0 ) .and. ( density <= 550.0) ) then

      arrhenius = tfm_density_arrhenius(1, ABDC, QBDC, temperature)

      strain_rate = strain_rate + (                                                &
      &  (-2.0_prec / 15.0_prec)                                                   &
      &  * (DELTA_B)                                                               &
      &  * ((8.0_prec * arrhenius * OMEGA) / (BOLTZMANN * temperature * (H**2.0))) &
      &  * (1.0_prec / (grain_radius * (MU**2.0)))                                 &
      &  * ((ICE_DENSITY / density)**3.0)                                          &
      &  * (1.0_prec - ((5.0_prec * density) / (3.0_prec * ICE_DENSITY)))          &
      &  * stress                                                                  &
      )
    
    else

      print *, 'module: tfm_density_processes'
      print *, 'function: tfm_density_grainBoundarySliding'
      print *, ''
      print *, 'The processes of grain boundary sliding is not defined '
      print *, 'for the given density:'
      print *, 'density = ', density
      print *, ''
      print *, 'Stopping right here!'
      STOP

    end if
  end function tfm_density_grainBoundarySliding


  function tfm_density_dislocationCreep(density, temperature, stress) &
    & result(strain_rate)
    implicit none

    real(prec), intent(in)    :: density
    real(prec), intent(in)    :: temperature
    real(prec), intent(in)    :: stress

    real(prec) :: strain_rate
    real(prec) :: rel_density
    real(prec) :: arrhenius
    real(prec) :: corrected_stress

    ! parameters
    real(prec) :: ADC = 3.22e-11_prec
    real(prec) :: QDC = 74500.0_prec

    arrhenius = tfm_density_arrhenius(1, ADC, QDC, temperature)
    rel_density = (density / ICE_DENSITY)

    if ( ( density > 550.0 ) .and. ( density <= 834.0 ) ) then
      
      strain_rate = (                                                           &
      &  (                                                                      &
      &    (-2.0_prec * arrhenius)                                              &
      &    * (1.0_prec - rel_density)                                           &
      &    * (((2.0_prec / ICE_N) * stress)**ICE_N)                             &
      &  )                                                                      &
      &  / ((1.0_prec - ((1.0_prec - rel_density)**(1.0_prec / ICE_N)))**ICE_N) &
      )

    else if ( ( density > 834.0 ) .and. ( density < ICE_DENSITY ) ) then
      
      corrected_stress = (                                            &
      &  stress - (AIR_PRESSURE * tfm_density_boyleMariotte(density)) &
      )
      
      strain_rate = (                                                             &
      &  ((-3.0_prec / 2.0_prec) * arrhenius)                                     &
      &  * (                                                                      &
      &    (1.0_prec - rel_density)                                               &
      &    / ((1.0_prec - ((1.0_prec - rel_density)**(1.0_prec / ICE_N)))**ICE_N) &
      &  )                                                                        &
      &  * (((3.0_prec / (2.0_prec * ICE_N)) * corrected_stress)**ICE_N)          &
      )

    else

      print *, 'module: tfm_density_processes'
      print *, 'function: tfm_density_dislocationCreep'
      print *, ''
      print *, 'The processes of dislocation creep is not defined for '
      print *, 'the given density:'
      print *, 'density = ', density
      print *, ''
      print *, 'Stopping right here!'
      STOP

    end if
  end function tfm_density_dislocationCreep


  function tfm_density_boundaryDiffusion(density, temperature, grain_radius, &
    & stress) result(strain_rate)
    implicit none

    real(prec), intent(in) :: density
    real(prec), intent(in) :: temperature
    real(prec), intent(in) :: grain_radius
    real(prec), intent(in) :: stress

    real(prec) :: strain_rate
    real(prec) :: arrhenius
    real(prec) :: corrected_stress

    ! parameters
    real(prec) :: ABDC = 3.0e-2_prec
    real(prec) :: QBDC = 44100.0_prec

    arrhenius = tfm_density_arrhenius(1, ABDC, QBDC, temperature)


    if ( ( density > 550.0 ) .and. ( density <= 834.0 ) ) then
      
      corrected_stress = stress
      
    else if ( ( density > 834.0 ) .and. ( density <= ICE_DENSITY ) ) then
      
      corrected_stress = (                                            &
      &  stress - (AIR_PRESSURE * tfm_density_boyleMariotte(density)) &
      )
      
    else

      print *, 'module: tfm_density_processes'
      print *, 'function: tfm_density_boundaryDiffusion'
      print *, ''
      print *, 'The processes of boundary diffusion is not defined for '
      print *, 'the given density:'
      print *, 'density = ', density
      print *, ''
      print *, 'Stopping right here!'
      STOP

    end if

    strain_rate = (                                             &
    &  (-37.0_prec / 2.0_prec)                                  &
    &  * (                                                      &
    &    (DELTA_B * arrhenius * OMEGA)                          &
    &    / (BOLTZMANN * temperature * (grain_radius**3.0_prec)) &
    &  )                                                        &
    &  * (ICE_DENSITY / density)                                &
    &  * corrected_stress                                       &
    )
  end function tfm_density_boundaryDiffusion


  function tfm_density_latticeDiffusion(density, temperature, grain_radius, &
    & stress) result(strain_rate)
    implicit none

    real(prec), intent(in) :: density
    real(prec), intent(in) :: temperature
    real(prec), intent(in) :: grain_radius
    real(prec), intent(in) :: stress
    
    real(prec) :: strain_rate
    real(prec) :: arrhenius
    real(prec) :: corrected_stress

    ! parameters
    real(prec), parameter :: ALD = 3.0e-2_prec
    real(prec), parameter :: QLD = 66200.0_prec


    arrhenius = tfm_density_arrhenius(1, ALD, QLD, temperature)

    if ( ( density > 550.0 ) .and. ( density <= 834.0 ) ) then
      
      strain_rate = (                                         &
      &  (-10.0_prec / 3.0_prec)                              &
      &  * (                                                  &
      &    (arrhenius * OMEGA)                                &
      &     / (BOLTZMANN * temperature * (grain_radius**2.0)) &
      &  )                                                    &
      &  * (ICE_DENSITY / density)                            &
      &  * stress                                             &
      )

    else if ( ( density > 834.0 ) .and. ( density <= ICE_DENSITY ) ) then
      
      corrected_stress = (                                            &
      &  stress - (AIR_PRESSURE * tfm_density_boyleMariotte(density)) &
      )
      
      strain_rate = (                                                                                &
      &  (-ICE_DENSITY / density)                                                                    &
      &  * ((3.0_prec * arrhenius * OMEGA) / (BOLTZMANN * temperature * grain_radius**2.0))          &
      &  * (grain_radius / (((ICE_DENSITY / (ICE_DENSITY - density))**(3.0_prec**-1.0)) - 1.0_prec)) &
      &  * ((ICE_DENSITY / density) * corrected_stress)                                              &
      )

    else

      print *, 'module: tfm_density_processes'
      print *, 'function: tfm_density_latticeDiffusion'
      print *, ''
      print *, 'The processes of lattice diffusion is not defined for '
      print *, 'the given density:'
      print *, 'density = ', density
      print *, ''
      print *, 'Stopping right here!'
      STOP
    
    end if
  end function tfm_density_latticeDiffusion
end module tfm_density_processes



module tfm_density_herronLangway
  use settings
  use tfm_constants
  implicit none

  contains


  subroutine tfm_density_mean_acc(nz, depth, density, age, mean_acc)
    implicit none

    integer, intent(in)                      :: nz
    real(prec), dimension(nz), intent(in)    :: depth
    real(prec), dimension(nz), intent(in)    :: density
    real(prec), dimension(nz), intent(in)    :: age
    real(prec), dimension(nz), intent(inout) :: mean_acc

    integer :: n

    mean_acc(nz) = 0.0
    do n = nz - 1, 1, -1
      mean_acc(n) = (depth(n+1) - depth(n)) * (density(n) / WATER_DENSITY)
      mean_acc(n) = mean_acc(n) + mean_acc(n+1)
    end do

    mean_acc(1:nz-1) = mean_acc(1:nz-1) / (age(1:nz-1) / SECONDS_YEAR)
  end subroutine tfm_density_mean_acc


  subroutine tfm_density_HLtype(nz, dt, stage1_params, &
    & stage2_params, depth, temperature, density, age, d_density)
    implicit none

    integer, intent(in)                   :: nz
    real(prec), intent(in)                :: dt
    real(prec), dimension(3), intent(in)  :: stage1_params
    real(prec), dimension(3), intent(in)  :: stage2_params
    real(prec), dimension(nz), intent(in) :: depth
    real(prec), dimension(nz), intent(in) :: temperature
    real(prec), dimension(nz), intent(in) :: density
    real(prec), dimension(nz), intent(in) :: age

    real(prec), dimension(nz), intent(inout) :: d_density

    real(prec), dimension(nz) :: mean_acc

    integer                   :: n
    real(prec), dimension(nz) :: c


    ! computation of the mean accumulation rate 
    ! over the lifetime of the firn parcel (kg a-1 m-2)
    call tfm_density_mean_acc(nz, depth, density, age, mean_acc)
    mean_acc = mean_acc * WATER_DENSITY

    ! boundary between first and seconds stage
    do n = nz, 1, -1
      
      ! first stage
      if ( ( density(n) > 0.0 ) .and. ( density(n) <= 550.0 ) ) then

        c(n) = (                                                   &
        &  stage1_params(1)                                        &
        &  * (mean_acc(n)**stage1_params(2))                       &
        &  * ACC_GRAVITY                                           &
        &  * exp(-stage1_params(3) / (GAS_CONST * temperature(n))) &
        )

      ! seconds stage
      else if ( ( density(n) > 550.0 ) .and. ( density(n) <= ICE_DENSITY ) ) then

        c(n) = (                                                   &
        &  stage2_params(1)                                        &
        &  * (mean_acc(n)**stage2_params(2))                       &
        &  * ACC_GRAVITY                                           &
        &  * exp(-stage2_params(3) / (GAS_CONST * temperature(n))) &
        )

      else

        print *, 'module: tfm_density, subroutine: tfm_density_HLtype'
        print *, 'The density seems to show an irregular value.'
        print *, 'Something went wrong! Stopping right here!'
        STOP
 
      end if
    end do
    
    ! density change
    d_density = ((dt / SECONDS_YEAR) * (c * (ICE_DENSITY - density)))
  end subroutine tfm_density_HLtype
end module tfm_density_herronLangway



module tfm_density_fischmeisterArzt
  use settings
  use tfm_constants
  implicit none

  interface tfm_density_arztCoordination
    module procedure tfm_density_arztCoordination_i, tfm_density_arztCoordination_n
  end interface tfm_density_arztCoordination

  interface tfm_density_arztContactarea
    module procedure tfm_density_arztContactarea_i, tfm_density_arztContactarea_n
  end interface tfm_density_arztContactarea

  ! parameters
  real(prec), parameter :: ZZERO = 7.3
  !real(prec), parameter :: ZZERO = 4.81
  real(prec), parameter :: CARZT = 15.5

  contains


  function tfm_density_arztCoordination_i(n, rel_density, d_zero) &
    & result(coordination_number)
    implicit none

    integer, intent(in)    :: n
    real(prec), intent(in) :: rel_density
    real(prec), intent(in) :: d_zero
    real(prec)             :: coordination_number
    integer                :: m

    m = 1 * n

    coordination_number = (                                          &
    &  ZZERO + CARZT * (((rel_density / d_zero)**(1.0 / 3.0)) - 1.0) &
    )
  end function tfm_density_arztCoordination_i


  function tfm_density_arztCoordination_n(n, rel_density, d_zero) &
    & result(coordination_number)
    implicit none

    integer, intent(in)                  :: n
    real(prec), dimension(n), intent(in) :: rel_density
    real(prec), intent(in)               :: d_zero
    real(prec), dimension(n)             :: coordination_number

    coordination_number = (                                          &
    &  ZZERO + CARZT * (((rel_density / d_zero)**(1.0 / 3.0)) - 1.0) &
    )
  end function tfm_density_arztCoordination_n


  function tfm_density_arztContactarea_i(n, rel_density, d_zero) &
    & result(contactarea)
    implicit none

    integer, intent(in)    :: n
    real(prec), intent(in) :: rel_density
    real(prec), intent(in) :: d_zero

    real(prec) :: coordination_number
    real(prec) :: r_i, r_ii
    real(prec) :: contactarea
    integer    :: m

    m = 1 * n

    coordination_number = tfm_density_arztCoordination(n, rel_density, d_zero)

    r_i = (rel_density / d_zero)**(1.0 / 3.0)
    r_ii = r_i + (                                                                    &
    &  (                                                                              &
    &    ((4.0 * ZZERO) * ((r_i - 1.0)**2.0) * ((2.0 * r_i) + 1.0))                   &
    &    + (CARZT * ((r_i - 1.0)**3.0) * ((3.0 * r_i) + 1.0))                         &
    &  )                                                                              &
    &  / (                                                                            &
    &    (12.0 * r_i)                                                                 &
    &    * ((4.0 * r_i) - (2.0 * ZZERO * (r_i - 1.0)) - (CARZT * ((r_i - 1.0)**2.0))) &
    &  )                                                                              &
    )

    contactarea = (                                     &
    &  (PI / (3.0 * coordination_number * (r_i**2.0)))  &
    &  * (                                              &
    &    (3.0 * ((r_ii**2.0) - 1.0) * ZZERO)            &
    &    + ((r_ii**2.0) * CARZT * ((2.0 * r_ii) - 3.0)) &
    &    + (CARZT)                                      &
    &  )                                                &
    )
  end function tfm_density_arztContactarea_i


  function tfm_density_arztContactarea_n(n, rel_density, d_zero) &
    & result(contactarea)
    implicit none

    integer, intent(in)                  :: n
    real(prec), dimension(n), intent(in) :: rel_density
    real(prec), intent(in)               :: d_zero

    real(prec), dimension(n) :: coordination_number
    real(prec), dimension(n) :: r_i, r_ii
    real(prec), dimension(n) :: contactarea

    coordination_number = tfm_density_arztCoordination(n, rel_density, d_zero)

    r_i = (rel_density / d_zero)**(1.0 / 3.0)
    r_ii = r_i + (                                                                    &
    &  (                                                                              &
    &    ((4.0 * ZZERO) * ((r_i - 1.0)**2.0) * ((2.0 * r_i) + 1.0))                   &
    &    + (CARZT * ((r_i - 1.0)**3.0) * ((3.0 * r_i) + 1.0))                         &
    &  )                                                                              &
    &  / (                                                                            &
    &    (12.0 * r_i)                                                                 &
    &    * ((4.0 * r_i) - (2.0 * ZZERO * (r_i - 1.0)) - (CARZT * ((r_i - 1.0)**2.0))) &
    &  )                                                                              &
    )

    contactarea = (                                     &
    &  (PI / (3.0 * coordination_number * (r_i**2.0)))  &
    &  * (                                              &
    &    (3.0 * ((r_ii**2.0) - 1.0) * ZZERO)            &
    &    + ((r_ii**2.0) * CARZT * ((2.0 * r_ii) - 3.0)) &
    &    + (CARZT)                                      &
    &  )                                                &
    )
  end function tfm_density_arztContactarea_n
end module tfm_density_fischmeisterArzt


module tfm_density_gagliardini
  use settings
  use tfm_constants
  use tfm_density_tools

  real(prec), parameter :: STAGE_DIV = 0.81_prec

  interface
    function invariant_inter(nz, param_a, param_b, strain_rate_inp) &
      & result(invariant)
      use settings
      implicit none

      integer, intent(in)                             :: nz
      real(prec), dimension(nz), intent(in)           :: param_a
      real(prec), dimension(nz), intent(in)           :: param_b
      real(prec), dimension(nz), intent(in), optional :: strain_rate_inp

      real(prec), dimension(nz) :: invariant
    end function invariant_inter
  end interface

  
  interface
    function viscosity_inter(nz, param, rate_factor, invariant) &
      & result(viscosity)
      use settings
      implicit none

      integer, intent(in)                   :: nz
      real(prec), dimension(nz), intent(in) :: param
      real(prec), dimension(nz), intent(in) :: rate_factor
      real(prec), dimension(nz), intent(in) :: invariant

      real(prec), dimension(nz) :: viscosity
    end function viscosity_inter
  end interface


  contains


  function tfm_density_gagliardiniParamA0(density) result(param_a0)
    implicit none

    real(prec), intent(in) :: density
    real(prec)             :: param_a0
    real(prec)             :: rel_density

    rel_density = density / ICE_DENSITY

    param_a0 = (                                         &
    &  (1.0 + ((2.0 / 3.0) * (1.0 - rel_density)))       &
    &  * (rel_density)**((-2.0 * ICE_N) / (ICE_N + 1.0)) &
    )
  end function tfm_density_gagliardiniParamA0


  function tfm_density_gagliardiniParamB0(density) result(param_b0)
    implicit none

    real(prec), intent(in) :: density
    real(prec)             :: param_b0
    real(prec)             :: rel_density
    
    rel_density = density / ICE_DENSITY

    param_b0 = (                                                  &
    &  (3.0 / 4.0)                                                &
    &  * (                                                        &
    &    ((1.0 - rel_density)**(1.0 / ICE_N))                     &
    &    / (ICE_N * (1.0 - ((1.0 - rel_density)**(1.0 / ICE_N)))) &
    &  )**((2.0 * ICE_N) / (ICE_N + 1.0))                         &
    )
  end function tfm_density_gagliardiniParamB0


  function tfm_density_gagliardiniRate(nz, temperature) &
    & result(rate_factor)
    implicit none

    integer, intent(in)                   :: nz
    real(prec), dimension(nz), intent(in) :: temperature

    integer                   :: n
    real(prec), dimension(nz) :: rate_factor

    ! parameters (values from Greve & Blatter, 2009)
    real(prec) :: TEMP_DIV = 263.15_prec
    real(prec) :: PRE_FACTOR_LOW  = 3.985e-13_prec
    real(prec) :: PRE_FACTOR_HIGH = 1.916e3_prec
    real(prec) :: ACTIVATION_ENERGY_LOW  =  60000.0_prec
    real(prec) :: ACTIVATION_ENERGY_HIGH = 139000.0_prec

    do n = 1, nz, 1
      if ( temperature(n) <= TEMP_DIV ) then
        
        rate_factor(n) = tfm_density_arrhenius( &
        &  1,                                   &
        &  PRE_FACTOR_LOW,                      &
        &  ACTIVATION_ENERGY_LOW,               &
        &  temperature(n)                       &
        &)
      else if ( temperature(n) > TEMP_DIV ) then

        rate_factor(n) = tfm_density_arrhenius( &
        &  1,                                   &
        &  PRE_FACTOR_HIGH,                     &
        &  ACTIVATION_ENERGY_HIGH,              &
        &  temperature(n)                       &
        &)

      else

        ! catch exception
        print *, 'module: tfm_density_gagliardini'
        print *, 'function: tfm_density_gagliardiniRate'
        print *, ''
        print *, 'It seems there are irregular temperature values!'
        print *, ''
        print *, 'Stopping right here!'
        STOP

      end if
    end do

    rate_factor = rate_factor**(-1.0_prec / ICE_N)
  end function tfm_density_gagliardiniRate

  
  function tfm_density_gagliardiniSolve(nz, density, stress, dt, param_a, &
    & param_b, rate_factor, invariant_func, shear_visco_func, bulk_visco_func) &
    & result(d_density)
    implicit none

    integer, intent(in)                   :: nz
    real(prec), intent(in)                :: dt
    real(prec), dimension(nz), intent(in) :: density
    real(prec), dimension(nz), intent(in) :: stress
    real(prec), dimension(nz), intent(in) :: param_a
    real(prec), dimension(nz), intent(in) :: param_b
    real(prec), dimension(nz), intent(in) :: rate_factor

    procedure(invariant_inter) :: invariant_func
    procedure(viscosity_inter) :: shear_visco_func
    procedure(viscosity_inter) :: bulk_visco_func

    real(prec), dimension(nz) :: d_density 

    integer                   :: n
    integer                   :: iter
    real(prec), dimension(nz) :: strain_rate
    real(prec), dimension(nz) :: d_density_prev
    real(prec), dimension(nz) :: invariant
    real(prec), dimension(nz) :: bulk_viscosity
    real(prec), dimension(nz) :: shear_viscosity

    integer, parameter :: MAX_ITER = 1000


    ! first guess of the invariant
    invariant = invariant_func(nz, param_a, param_b)

    ! viscosities
    shear_viscosity = shear_visco_func(    &
    &  nz, param_a, rate_factor, invariant &
    &)
    bulk_viscosity = bulk_visco_func(      &
    &  nz, param_b, rate_factor, invariant & 
    &)

    ! strain rate
    strain_rate = (                                                             &
    &  (1.0_prec/ (((4.0_prec / 3.0_prec) * shear_viscosity) + bulk_viscosity)) &
    &  * stress                                                                 &
    &)

    d_density      = -999999.9
    d_density_prev = +999999.9
    iter           = 0

    do while ( (maxval(abs(d_density - d_density_prev)) > 1.0e-2) .or. (iter == MAX_ITER) )

      ! first guess of the invariant
      invariant = invariant_func( &
      &  nz, param_a, param_b, strain_rate          &
      &)

      ! viscosities
      shear_viscosity = shear_visco_func(    &
      &  nz, param_a, rate_factor, invariant &
      &)
      bulk_viscosity = bulk_visco_func(      &
      &  nz, param_b, rate_factor, invariant &
      &)

      ! strain rate
      strain_rate = (                                                             &
      &  (1.0_prec/ (((4.0_prec / 3.0_prec) * shear_viscosity) + bulk_viscosity)) &
      &  * stress                                                                 &
      &)

      ! densification
      d_density_prev = 1.0_prec * d_density
      d_density = dt * strain_rate * density

      iter = iter + 1
    end do

    ! avoid the singularity at ice density
    do n = 1, nz, 1
      if ( (density(n) + d_density(n)) > (ICE_DENSITY - 10.0e-5) ) then
        d_density(n) = (ICE_DENSITY - density(n)) - 10.0e-5
      end if
    end do
  end function tfm_density_gagliardiniSolve
end module tfm_density_gagliardini



module tfm_density
  use settings
  use tfm_constants
  use tfm_density_tools
  use tfm_density_herronLangway
  use tfm_density_stress
  use tfm_density_fischmeisterArzt
  use tfm_density_processes
  use tfm_density_gagliardini
  implicit none


  contains


  function tfm_density_depth(nz, depth, density, d_density) result(d_depth)
    implicit none

    integer, intent(in)                   :: nz
    real(prec), dimension(nz), intent(in) :: depth
    real(prec), dimension(nz), intent(in) :: density
    real(prec), dimension(nz), intent(in) :: d_density

    real(prec), dimension(nz) :: d_depth

    real(prec), dimension(nz) :: dz
    real(prec), dimension(nz) :: ddz
    integer                   :: n

    dz(1) = 0.0_prec
    dz(2:nz) = (depth(2:nz) - depth(1:nz-1))
    ddz = (dz * (density / (density + d_density))) - dz

    d_depth(1) = 0.0_prec
    do n = 2, nz, 1
      d_depth(n) = d_depth(n-1) + ddz(n)
    end do
  end function tfm_density_depth


  function tfm_density_gagliardini1998(nz, dt, depth, density, temperature, &
    & age, grain_radius) result(d_density)
    implicit none

    integer, intent(in)                   :: nz
    real(prec), intent(in)                :: dt
    real(prec), dimension(nz), intent(in) :: depth
    real(prec), dimension(nz), intent(in) :: density
    real(prec), dimension(nz), intent(in) :: temperature
    real(prec), dimension(nz), intent(in) :: age
    real(prec), dimension(nz), intent(in) :: grain_radius

    integer                   :: n
    real(prec)                :: rel_density
    real(prec), dimension(nz) :: d_density
    real(prec), dimension(nz) :: param_a
    real(prec), dimension(nz) :: param_b
    real(prec), dimension(nz) :: rate_factor
    real(prec), dimension(nz) :: stress

    ! doing nothing
    call tfm_density_do_nothing(nz, age)
    call tfm_density_do_nothing(nz, grain_radius)

    ! stress
    call tfm_density_computeStress(nz, depth, density, stress)

    ! density dependent parameter of the model as defined by Greve & BLatter 2009
    do n = 1, nz, 1
      
      rel_density = density(n) / ICE_DENSITY
      
      !if ( rel_density < 0.5 ) then
      !  param_b(n) = exp(                &
      !  &  (451.63 * (rel_density**2.0)) &
      !  &  - (474.34 * rel_density)      &
      !  &  + 128.12                      &
      !  )

      !else if ( (rel_density >= 0.5) .and. (rel_density < 0.785) ) then
      !  param_b(n) = exp((-17.15 * rel_density) + 12.42)

      !else if ( (rel_density >= 0.785) .and. (rel_density < 1.0) ) then
      !  param_b(n) = tfm_density_gagliardiniParamB0(density(n))

      if ( (rel_density > 0.0) .and. (rel_density < 0.785) ) then
        param_a(n) = exp((-19.67 * rel_density) + 15.94)
        param_b(n) = exp((-27.65 * rel_density) + 20.37)

      else if ( (rel_density >= 0.785) .and. (rel_density < 1.0) ) then
        param_a(n) = tfm_density_gagliardiniParamA0(density(n))
        param_b(n) = tfm_density_gagliardiniParamB0(density(n))

      else
        ! catch exception
        print *, 'module: tfm_density'
        print *, 'function: tfm_density_gagliardini1998'
        print *, ''
        print *, 'It seems the density exceeds the range of valid '
        print *, 'values at some point!'
        print *, ''
        print *, 'Stopping right here!'
        STOP

      end if

      !param_a(n) = (                                    &
      !&  param_b(n) * (                                 &
      !&    tfm_density_gagliardiniParamA0(density(n))   &
      !&    / tfm_density_gagliardiniParamB0(density(n)) &
      !&  )                                              &
      !)

    end do

    ! temperature dependent associated rate factor
    rate_factor = tfm_density_gagliardiniRate(nz, temperature)

    ! solving for the density change
    d_density = tfm_density_gagliardiniSolve(                  &
    &  nz, density, stress, dt, param_a, param_b, rate_factor, &
    &  invariant_func, visco_func, visco_func                  &
    )

    contains

    function invariant_func(nz, param_a, param_b, &
      & strain_rate_inp) result(invariant)
      implicit none

      integer, intent(in)                             :: nz
      real(prec), dimension(nz), intent(in)           :: param_a
      real(prec), dimension(nz), intent(in)           :: param_b
      real(prec), dimension(nz), intent(in), optional :: strain_rate_inp

      real(prec), dimension(nz) :: invariant
      real(prec), dimension(nz) :: strain_rate

      if ( present(strain_rate_inp) ) then
        strain_rate = strain_rate_inp
      else
        strain_rate = 1.0e-10_prec
      end if

      invariant = (                                           &
      &  strain_rate                                          & 
      &  * (((3.0 / (4.0 * param_a)) + (1.0 / param_b))**0.5) &
      )
    end function invariant_func


    function visco_func(nz, param, rate_factor, invariant) &
      & result(viscosity)

      integer, intent(in) :: nz
      real(prec), dimension(nz), intent(in) :: param
      real(prec), dimension(nz), intent(in) :: rate_factor
      real(prec), dimension(nz), intent(in) :: invariant

      real(prec), dimension(nz) :: viscosity

      viscosity = (                                        &
      &  (1.0_prec / param)                                &
      &  * rate_factor                                     &
      &  * (invariant**(-(1.0_prec - (1.0_prec / ICE_N)))) &
      )
    end function visco_func
  end function tfm_density_gagliardini1998


  function tfm_density_timmsfit(nz, dt, depth, density, temperature, &
    & age, grain_radius) result(d_density)
    implicit none

    integer, intent(in)                   :: nz
    real(prec), intent(in)                :: dt
    real(prec), dimension(nz), intent(in) :: depth
    real(prec), dimension(nz), intent(in) :: density
    real(prec), dimension(nz), intent(in) :: temperature
    real(prec), dimension(nz), intent(in) :: age
    real(prec), dimension(nz), intent(in) :: grain_radius

    integer                   :: n
    real(prec)                :: rel_density
    real(prec), dimension(nz) :: d_density
    real(prec), dimension(nz) :: param_a
    real(prec), dimension(nz) :: param_b
    real(prec), dimension(nz) :: rate_factor
    real(prec), dimension(nz) :: stress


    ! doing nothing
    call tfm_density_do_nothing(nz, age)
    call tfm_density_do_nothing(nz, grain_radius)

    ! stress
    call tfm_density_computeStress(nz, depth, density, stress)

    ! density dependent parameter of the model as defined by Greve & BLatter 2009
    do n = 1, nz, 1
      rel_density = density(n) / ICE_DENSITY

      if ( (rel_density > 0.0) .and. (rel_density <= 0.79) ) then
        param_a(n) = exp(                 &
        &  24.60215                       &
        &  - (58.573530 * rel_density)    &
        &  - (-35.5 * (rel_density**2.0)) &
        )
        param_b(n) = (                                    &
        &  (                                              &
        &    tfm_density_gagliardiniParamB0(density(n))   &
        &    / tfm_density_gagliardiniParamA0(density(n)) &
        &  ) * param_a(n)                                 &
        )

      else if ( (rel_density > 0.79) .and. (rel_density < 1.0) ) then
        param_a(n) = tfm_density_gagliardiniParamA0(density(n))
        param_b(n) = tfm_density_gagliardiniParamB0(density(n))

      else
        ! catch exception
        print *, 'module: tfm_density'
        print *, 'function: tfm_density_timmsfit'
        print *, ''
        print *, 'It seems the density exceeds the range of valid '
        print *, 'values at some point!'
        print *, ''
        print *, 'Stopping right here!'
        STOP

      end if
    end do

    ! temperature dependent associated rate factor
    rate_factor = tfm_density_gagliardiniRate(nz, temperature)

    ! solving for the density change
    d_density = tfm_density_gagliardiniSolve(                  &
    &  nz, density, stress, dt, param_a, param_b, rate_factor, &
    &  invariant_func, visco_func, visco_func                  &
    )

    contains

    function invariant_func(nz, param_a, param_b, &
      & strain_rate_inp) result(invariant)
      implicit none

      integer, intent(in)                             :: nz
      real(prec), dimension(nz), intent(in)           :: param_a
      real(prec), dimension(nz), intent(in)           :: param_b
      real(prec), dimension(nz), intent(in), optional :: strain_rate_inp

      real(prec), dimension(nz) :: invariant
      real(prec), dimension(nz) :: strain_rate

      if ( present(strain_rate_inp) ) then
        strain_rate = strain_rate_inp
      else
        strain_rate = 1.0e-10_prec
      end if

      invariant = (                                                   &
      &  strain_rate                                                  &
      &  * (((1.0 / (3.0 * param_a)) + (1.0 / (4.0 * param_b)))**0.5) &
      )
    end function invariant_func


    function visco_func(nz, param, rate_factor, invariant) &
      & result(viscosity)

      integer, intent(in) :: nz
      real(prec), dimension(nz), intent(in) :: param
      real(prec), dimension(nz), intent(in) :: rate_factor
      real(prec), dimension(nz), intent(in) :: invariant

      real(prec), dimension(nz) :: viscosity

      viscosity = (                                        &
      &  (1.0_prec / (2.0 * param))                        &
      &  * rate_factor                                     &
      &  * (invariant**(-(1.0_prec - (1.0_prec / ICE_N)))) &
      )
    end function visco_func
  end function tfm_density_timmsfit


  function tfm_density_greve2009(nz, dt, depth, density, temperature, &
    & age, grain_radius) result(d_density)
    implicit none

    integer, intent(in)                   :: nz
    real(prec), intent(in)                :: dt
    real(prec), dimension(nz), intent(in) :: depth
    real(prec), dimension(nz), intent(in) :: density
    real(prec), dimension(nz), intent(in) :: temperature
    real(prec), dimension(nz), intent(in) :: age
    real(prec), dimension(nz), intent(in) :: grain_radius

    integer                   :: n
    real(prec)                :: rel_density
    real(prec), dimension(nz) :: d_density
    real(prec), dimension(nz) :: param_a
    real(prec), dimension(nz) :: param_b
    real(prec), dimension(nz) :: rate_factor
    real(prec), dimension(nz) :: stress


    ! doing nothing
    call tfm_density_do_nothing(nz, age)
    call tfm_density_do_nothing(nz, grain_radius)

    ! stress
    call tfm_density_computeStress(nz, depth, density, stress)

    ! density dependent parameter of the model as defined by Greve & BLatter 2009
    do n = 1, nz, 1
      rel_density = density(n) / ICE_DENSITY

      if ( (rel_density > 0.0) .and. (rel_density <= 0.81) ) then
        param_a(n) = exp(13.22240 - (15.78652 * rel_density))
        param_b(n) = exp(15.09371 - (20.46489 * rel_density))

      else if ( (rel_density > 0.81) .and. (rel_density < 1.0) ) then
        param_a(n) = tfm_density_gagliardiniParamA0(density(n))
        param_b(n) = tfm_density_gagliardiniParamB0(density(n))

      else
        ! catch exception
        print *, 'module: tfm_density'
        print *, 'function: tfm_density_greve2009'
        print *, ''
        print *, 'It seems the density exceeds the range of valid '
        print *, 'values at some point!'
        print *, ''
        print *, 'Stopping right here!'
        STOP

      end if
    end do

    ! temperature dependent associated rate factor
    rate_factor = tfm_density_gagliardiniRate(nz, temperature)

    ! solving for the density change
    d_density = tfm_density_gagliardiniSolve(                  &
    &  nz, density, stress, dt, param_a, param_b, rate_factor, &
    &  invariant_func, visco_func, visco_func                  &
    )

    contains

    function invariant_func(nz, param_a, param_b, &
      & strain_rate_inp) result(invariant)
      implicit none

      integer, intent(in)                             :: nz
      real(prec), dimension(nz), intent(in)           :: param_a
      real(prec), dimension(nz), intent(in)           :: param_b
      real(prec), dimension(nz), intent(in), optional :: strain_rate_inp

      real(prec), dimension(nz) :: invariant
      real(prec), dimension(nz) :: strain_rate

      if ( present(strain_rate_inp) ) then
        strain_rate = strain_rate_inp
      else
        strain_rate = 1.0e-10_prec
      end if

      invariant = (                                                   &
      &  strain_rate                                                  &
      &  * (((1.0 / (3.0 * param_a)) + (1.0 / (4.0 * param_b)))**0.5) &
      )
    end function invariant_func


    function visco_func(nz, param, rate_factor, invariant) &
      & result(viscosity)

      integer, intent(in) :: nz
      real(prec), dimension(nz), intent(in) :: param
      real(prec), dimension(nz), intent(in) :: rate_factor
      real(prec), dimension(nz), intent(in) :: invariant

      real(prec), dimension(nz) :: viscosity

      viscosity = (                                        &
      &  (1.0_prec / (2.0 * param))                        &
      &  * rate_factor                                     &
      &  * (invariant**(-(1.0_prec - (1.0_prec / ICE_N)))) &
      )
    end function visco_func
  end function tfm_density_greve2009


  function tfm_density_zwinger2007(nz, dt, depth, density, temperature, &
    & age, grain_radius) result(d_density)
    implicit none

    integer, intent(in)                   :: nz
    real(prec), intent(in)                :: dt
    real(prec), dimension(nz), intent(in) :: depth
    real(prec), dimension(nz), intent(in) :: density
    real(prec), dimension(nz), intent(in) :: temperature
    real(prec), dimension(nz), intent(in) :: age
    real(prec), dimension(nz), intent(in) :: grain_radius

    integer                   :: n
    real(prec)                :: rel_density
    real(prec), dimension(nz) :: d_density
    real(prec), dimension(nz) :: param_a
    real(prec), dimension(nz) :: param_b
    real(prec), dimension(nz) :: rate_factor
    real(prec), dimension(nz) :: stress

    ! doing nothing
    call tfm_density_do_nothing(nz, age)
    call tfm_density_do_nothing(nz, grain_radius)

    ! stress
    call tfm_density_computeStress(nz, depth, density, stress)

    ! density dependent parameter of the model as defined by Greve & BLatter 2009
    do n = 1, nz, 1
      rel_density = density(n) / ICE_DENSITY

      if ( (rel_density > 0.0) .and. (rel_density <= 0.81) ) then
        param_a(n) = exp(13.22240 - (15.78652 * rel_density))
        param_b(n) = exp(15.09371 - (20.46489 * rel_density))

      else if ( (rel_density > 0.81) .and. (rel_density < 1.0) ) then
        param_a(n) = tfm_density_gagliardiniParamA0(density(n))
        param_b(n) = tfm_density_gagliardiniParamB0(density(n))
 
      else
        ! catch exception
        print *, 'module: tfm_density'
        print *, 'function: tfm_density_zwinger2007'
        print *, ''
        print *, 'It seems the density exceeds the range of valid '
        print *, 'values at some point!'
        print *, ''
        print *, 'Stopping right here!'
        STOP

      end if
    end do

    ! temperature dependent associated rate factor
    rate_factor = tfm_density_gagliardiniRate(nz, temperature)

    ! solving for the density change
    d_density = tfm_density_gagliardiniSolve(                  &
    &  nz, density, stress, dt, param_a, param_b, rate_factor, &
    &  invariant_func, shear_visco_func, bulk_visco_func       &
    )

    contains

    function invariant_func(nz, param_a, param_b, &
      & strain_rate_inp) result(invariant)
      implicit none

      integer, intent(in)                             :: nz
      real(prec), dimension(nz), intent(in)           :: param_a
      real(prec), dimension(nz), intent(in)           :: param_b
      real(prec), dimension(nz), intent(in), optional :: strain_rate_inp

      real(prec), dimension(nz) :: invariant
      real(prec), dimension(nz) :: strain_rate

      if ( present(strain_rate_inp) ) then
        strain_rate = strain_rate_inp
      else
        strain_rate = 1.0e-10_prec
      end if

      invariant = (                                           &
      &  strain_rate                                          &
      &  * (((3.0 / (4.0 * param_a)) + (1.0 / param_b))**0.5) &
      )
    end function invariant_func


    function shear_visco_func(nz, param_a, rate_factor, invariant) &
      & result(shear_viscosity)

      integer, intent(in) :: nz
      real(prec), dimension(nz), intent(in) :: param_a
      real(prec), dimension(nz), intent(in) :: rate_factor
      real(prec), dimension(nz), intent(in) :: invariant

      real(prec), dimension(nz) :: shear_viscosity

      shear_viscosity = (                                  &
      &  (2.0_prec / param_a)                              &
      &  * rate_factor                                     &
      &  * (invariant**(-(1.0_prec - (1.0_prec / ICE_N)))) &
      )
    end function shear_visco_func


    function bulk_visco_func(nz, param_b, rate_factor, invariant) &
      & result(bulk_viscosity)

      integer, intent(in) :: nz
      real(prec), dimension(nz), intent(in) :: param_b
      real(prec), dimension(nz), intent(in) :: rate_factor
      real(prec), dimension(nz), intent(in) :: invariant

      real(prec), dimension(nz) :: bulk_viscosity

      bulk_viscosity = (                                   &
      &  (1.0_prec / param_b)                              &
      &  * rate_factor                                     &
      &  * (invariant**(-(1.0_prec - (1.0_prec / ICE_N)))) &
      )
    end function bulk_visco_func
  end function tfm_density_zwinger2007


  function tfm_density_breant2017(nz, dt, depth, density, temperature, &
    & age, grain_radius) result(d_density)
    implicit none

    integer, intent(in)                   :: nz
    real(prec), intent(in)                :: dt
    real(prec), dimension(nz), intent(in) :: depth
    real(prec), dimension(nz), intent(in) :: density
    real(prec), dimension(nz), intent(in) :: temperature
    real(prec), dimension(nz), intent(in) :: age
    real(prec), dimension(nz), intent(in) :: grain_radius

    integer                   :: n, m
    integer                   :: mz
    real(prec), dimension(nz) :: rel_density
    real(prec), dimension(nz) :: d_density
    real(prec), dimension(nz) :: arrhenius_dc
    real(prec), dimension(nz) :: stress
    real(prec), dimension(nz) :: eff_stress
    real(prec), dimension(nz) :: coord_number
    real(prec), dimension(nz) :: contact_area
    real(prec), dimension(nz) :: strain_rate
    real(prec)                :: gamma_gbs 
    real(prec)                :: dz

    real(prec) :: arrhenius_gamma
    real(prec) :: contact_area_gamma
    real(prec) :: eff_stress_gamma 
    real(prec) :: stress_gamma
    real(prec) :: temperature_gamma

    real(prec), dimension(:,:), allocatable :: bagmean
    
    real(prec), parameter :: DZERO = 0.56
    real(prec), parameter :: MIN_STRESS = 0.1e5
    real(prec), parameter :: A0 = 7.89e-15
    real(prec), parameter :: A1 = 1.05e9,  Q1 = 110000.0
    real(prec), parameter :: A2 = 1400.0,  Q2 =  75000.0
    real(prec), parameter :: A3 = 6.0e-15, Q3 =   1500.0
    real(prec), parameter :: QGBS = 49500.0

    call tfm_density_do_nothing(nz, age)
    call tfm_density_do_nothing(nz, grain_radius)

    rel_density = density / ICE_DENSITY

    ! temperature dependence
    arrhenius_dc = A0 * (                            &
    &    (A1 * exp(-Q1 / (GAS_CONST * temperature))) &
    &  + (A2 * exp(-Q2 / (GAS_CONST * temperature))) &
    &  + (A3 * exp(-Q3 / (GAS_CONST * temperature))) &
    )

    ! model by Arzt
    coord_number = tfm_density_arztCoordination(nz, rel_density, DZERO)
    contact_area = tfm_density_arztContactarea(nz, rel_density, DZERO)

    ! stress
    call tfm_density_computeStress(nz, depth, density, stress)
    do n = 1, nz, 1
      if (stress(n) < MIN_STRESS) stress(n) = MIN_STRESS
    end do
    eff_stress = stress * (                                     &
    &  (4.0 * PI) / (contact_area * coord_number * rel_density) &
    )


    dz = 1.0
    mz = floor((depth(nz) - depth(1)) / dz)
    allocate(bagmean(2,mz))
    call tfm_density_bagmean(nz, depth, rel_density, dz, mz, bagmean)

    do m = mz, 1, -1
      if ( bagmean(2,m) > 0.6 ) EXIT
    end do

    do n = nz, 1, -1
      if ( depth(n) < bagmean(1,m+1) ) EXIT
    end do

    call tfm_density_lin_interp(      &
    &  depth(n+1), depth(n),          &
    &  stress(n+1), stress(n),        &
    &  (bagmean(1,m+1) - depth(n+1)), &
    &  stress_gamma                   &
    )
    call tfm_density_lin_interp(         &
    &  depth(n+1), depth(n),             &
    &  temperature(n+1), temperature(n), &
    &  (bagmean(1,m+1) - depth(n+1)),    &
    &  temperature_gamma                 &
    )

    do n = nz, 1, -1
      if ( depth(n) < bagmean(1,m) ) EXIT
    end do

    call tfm_density_lin_interp(           &
    &  depth(n+1), depth(n),               &
    &  arrhenius_dc(n+1), arrhenius_dc(n), &
    &  (bagmean(1,m) - depth(n+1)),        &
    &  arrhenius_gamma                     &
    )
    call tfm_density_lin_interp(           &
    &  depth(n+1), depth(n),               &
    &  contact_area(n+1), contact_area(n), &
    &  (bagmean(1,m) - depth(n+1)),        &
    &  contact_area_gamma                  &
    )
    call tfm_density_lin_interp(           &
    &  depth(n+1), depth(n),               &
    &  eff_stress(n+1), eff_stress(n),     &
    &  (bagmean(1,m) - depth(n+1)),        &
    &  eff_stress_gamma                    &
    )

    gamma_gbs = (                                     &
    &  (5.3 * arrhenius_gamma)                        &
    &  * ((DZERO * (bagmean(2,m)**2.0))**(1.0 / 3.0)) &
    &  * ((contact_area_gamma / PI)**0.5)             &
    &  * ((eff_stress_gamma / 3.0)**ICE_N)            &
    )
    gamma_gbs = gamma_gbs / (                                 &
    &  (stress_gamma / (bagmean(2,m+1)**2.0))                 &
    &  * (1.0 + (0.5 / 6.0) - ((5.0 / 3.0) * bagmean(2,m+1))) &
    &  * (exp(-QGBS / (GAS_CONST * temperature_gamma)))       &
    )


    do n = 1, nz, 1

      if ( rel_density(n) < 0.6 ) then

        strain_rate(n) = (                                        &
        &  gamma_gbs                                              &
        &  * (stress(n) / (rel_density(n)**2.0))                  &
        &  * (1.0 + (0.5 / 6.0) - ((5.0 / 3.0) * rel_density(n))) &
        &  * (exp(-QGBS / (GAS_CONST * temperature(n))))          &
        )

      else if ( ( rel_density(n) >= 0.6 ) .and. ( rel_density(n) < 1.0) ) then

        if ( rel_density(n) > 0.9 ) then
          eff_stress(n) = tfm_density_rel_boylemariotte(rel_density(n))
        end if

        strain_rate(n) = (                                  &
        &  (5.3 * arrhenius_dc(n))                          &
        &  * ((DZERO * (rel_density(n)**2.0))**(1.0 / 3.0)) &
        &  * ((contact_area(n) / PI)**0.5)                  &
        &  * ((eff_stress(n) / 3.0)**ICE_N)                 &
        )

      else

        ! catch exception
        print *, 'module: tfm_density, function tfm_density_brean2017'
        print *, 'The density exceeds ice density!'
        print *, 'Stopping right here!'
        STOP

      end if
    end do

    ! densification
    d_density = dt * strain_rate * density
  end function tfm_density_breant2017


  function tfm_density_medley2020(nz, dt, depth, density, temperature, &
    & age, grain_radius) result(d_density)
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
    real(prec), dimension(nz), intent(in) :: depth
    real(prec), dimension(nz), intent(in) :: density
    real(prec), dimension(nz), intent(in) :: temperature
    real(prec), dimension(nz), intent(in) :: age
    real(prec), dimension(nz), intent(in) :: grain_radius

    real(prec), dimension(nz) :: d_density

    integer                   :: n
    real(prec)                :: mean_temperature

    call tfm_density_do_nothing(nz, grain_radius)

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
    &  (/ a0, ALPHA0, EC0 /), &
    &  (/ a1, ALPHA1, EC1 /), &
    &  depth,                 &
    &  temperature,           &
    &  density,               &
    &  age,                   &
    &  d_density              &
    )
  end function tfm_density_medley2020 


  function tfm_density_herron1980(nz, dt, depth, density, temperature, &
    & age, grain_radius) result(d_density)
    implicit none

    ! first stage parameters
    real(prec), parameter :: ALPHA0 = 1.0
    real(prec), parameter :: EC0    = 10160.0
    real(prec), parameter :: A0     = 11.0 * (0.001**ALPHA0) / ACC_GRAVITY

    ! second stage parameters
    real(prec), parameter :: ALPHA1 = 0.5
    real(prec), parameter :: EC1    = 21400.0
    real(prec), parameter :: A1     = 575.0 * (0.001**ALPHA1) / ACC_GRAVITY

    integer, intent(in)                   :: nz
    real(prec), intent(in)                :: dt
    real(prec), dimension(nz), intent(in) :: depth
    real(prec), dimension(nz), intent(in) :: density
    real(prec), dimension(nz), intent(in) :: temperature
    real(prec), dimension(nz), intent(in) :: age
    real(prec), dimension(nz), intent(in) :: grain_radius

    real(prec), dimension(nz) :: d_density

    call tfm_density_do_nothing(nz, grain_radius)

    ! call Herron & Langway model
    call tfm_density_HLtype(  &
    &  nz, dt,                &
    &  (/ A0, ALPHA0, EC0 /), &
    &  (/ A1, ALPHA1, EC1 /), &
    &  depth,                 &
    &  temperature,           &
    &  density,               &
    &  age,                   &
    &  d_density              &
    )
  end function tfm_density_herron1980


  function tfm_density_arthern1998(nz, dt, depth, density, temperature, &
    & age, grain_radius) result(d_density)
    implicit none

    integer, intent(in)                   :: nz
    real(prec), intent(in)                :: dt
    real(prec), dimension(nz), intent(in) :: depth
    real(prec), dimension(nz), intent(in) :: density
    real(prec), dimension(nz), intent(in) :: temperature
    real(prec), dimension(nz), intent(in) :: age
    real(prec), dimension(nz), intent(in) :: grain_radius

    integer                   :: n
    real(prec), dimension(nz) :: d_density
    real(prec), dimension(nz) :: stress
    real(prec), dimension(nz) :: strain_rate

    call tfm_density_do_nothing(nz, age)

    ! computation of the current stress state
    call tfm_density_computeStress(nz, depth, density, stress)

    ! compuation of the strain rate
    strain_rate(:) = 0.0_prec
    do n = 1, nz, 1
      if ( ( density(n) > 0.0 ) .and. ( density(n) <= 550.0 ) ) then

        strain_rate(n) = strain_rate(n) + (  &
        &  tfm_density_grainBoundarySliding( &
        &    density(n),                     &
        &    temperature(n),                 &
        &    grain_radius(n),                &
        &    stress(n)                       &
        &  )                                 &
        )

      else if ( ( density(n) > 550.0 ) .and. ( density(n) < 917.0 ) ) then
        
        strain_rate(n) = strain_rate(n) + ( &
        &  tfm_density_dislocationCreep(    &
        &    density(n),                    &
        &    temperature(n),                &
        &    stress(n)                      &
        &  )                                &
        )

        strain_rate(n) = strain_rate(n) + ( &
        &  tfm_density_boundaryDiffusion(   &
        &    density(n),                    &
        &    temperature(n),                &
        &    grain_radius(n),               &
        &    stress(n)                      &
        &  )                                &
        )

        strain_rate(n) = strain_rate(n) + ( &
        &  tfm_density_latticeDiffusion(    &
        &    density(n),                    &
        &    temperature(n),                &
        &    grain_radius(n),               &
        &    stress(n)                      &
        &  )                                &
        )

      else if ( density(n) >= 917.0 ) then

         strain_rate(n) = strain_rate(n) + 0.0_prec

      else
        print *, 'Module: tfm_density, Function: tfm_density_arthern1998'
        print *, 'The density seems to show an irregular value.'
        print *, 'Something went wrong! Stopping right here!'
        STOP
      end if
    end do

    ! densification from strain rate
    d_density = -dt * strain_rate * density
  end function tfm_density_arthern1998


  function tfm_density_li2003(nz, dt, depth, density, temperature, &
    & age, grain_radius) result(d_density)
    implicit none

    ! first stage parameters
    real(prec), parameter :: ALPHA0 = 1.0
    real(prec)            :: ec0
    real(prec)            :: a0

    ! second stage parameters
    real(prec), parameter :: ALPHA1 = 1.0
    real(prec)            :: ec1
    real(prec)            :: a1

    integer, intent(in)                   :: nz
    real(prec), intent(in)                :: dt
    real(prec), dimension(nz), intent(in) :: depth
    real(prec), dimension(nz), intent(in) :: density
    real(prec), dimension(nz), intent(in) :: temperature
    real(prec), dimension(nz), intent(in) :: age
    real(prec), dimension(nz), intent(in) :: grain_radius

    integer                   :: n
    real(prec), dimension(nz) :: d_density
    real(prec)                :: mean_temperature
    real(prec)                :: grain_growth_rate
    real(prec)                :: beta

    call tfm_density_do_nothing(nz, grain_radius)

    ! 10 m temperature
    do n = nz, 1, -1
      if ( (depth(n) - depth(nz)) <= -10.0 ) EXIT
    end do
    mean_temperature = temperature(n)

    ! grain growth rate and beta factor
    grain_growth_rate = (                                     &
    &  8.36_prec                                              &
    &  * (abs(mean_temperature - TEMP_OFFSET)**(-2.061_prec)) &
    )
    beta = 139.21_prec - (0.542_prec * mean_temperature)

    ! translation to a HL-model
    a0 = (grain_growth_rate * beta) / ICE_DENSITY
    a1 = a0

    ! activation energies
    ec0 = (                                                   &
    &  883.8_prec                               &
    &  * (abs(mean_temperature - TEMP_OFFSET)**(-0.885_prec)) &
    )
    ec1 = ec0

    ! call Herron & Langway model
    call tfm_density_HLtype(  &
    &  nz, dt,                &
    &  (/ a0, ALPHA0, ec0 /), &
    &  (/ a1, ALPHA1, ec1 /), &
    &  depth,                 &
    &  temperature,           &
    &  density,               &
    &  age,                   &
    &  d_density              &
    )
  end function tfm_density_li2003


  function tfm_density_helsen2008(nz, dt, depth, density, temperature, &
    & age, grain_radius) result(d_density)
    implicit none

    ! first stage parameters
    real(prec), parameter :: ALPHA0 = 1.0
    real(prec), parameter :: EC0    = 0.0
    real(prec)            :: a0

    ! second stage parameters
    real(prec), parameter :: ALPHA1 = 1.0
    real(prec), parameter :: EC1    = 0.0
    real(prec)            :: a1

    integer, intent(in)                   :: nz
    real(prec), intent(in)                :: dt
    real(prec), dimension(nz), intent(in) :: depth
    real(prec), dimension(nz), intent(in) :: density
    real(prec), dimension(nz), intent(in) :: temperature
    real(prec), dimension(nz), intent(in) :: age
    real(prec), dimension(nz), intent(in) :: grain_radius

    integer                   :: n
    real(prec), dimension(nz) :: d_density
    real(prec)                :: mean_temperature
    real(prec)                :: grain_growth_rate
    real(prec)                :: beta

    call tfm_density_do_nothing(nz, grain_radius)

    ! 10 m temperature
    do n = nz, 1, -1
      if ( (depth(n) - depth(nz)) <= -10.0 ) EXIT
    end do
    mean_temperature = temperature(n)

    ! grain growth rate and beta factor
    grain_growth_rate = (                                     &
    &  8.36_prec                                              &
    &  * (abs(mean_temperature - TEMP_OFFSET)**(-2.061_prec)) &
    )
    beta = 76.138_prec - (0.28965_prec * mean_temperature)

    ! translation to a HL-model
    a0 = (grain_growth_rate * beta) / ICE_DENSITY
    a1 = (grain_growth_rate * beta) / ICE_DENSITY

    ! call Herron & Langway model
    call tfm_density_HLtype(  &
    &  nz, dt,                &
    &  (/ a0, ALPHA0, EC0 /), &
    &  (/ a1, ALPHA1, EC1 /), &
    &  depth,                 &
    &  temperature,           &
    &  density,               &
    &  age,                   &
    &  d_density              &
    )
  end function tfm_density_helsen2008


  function tfm_density_arthern2010(nz, dt, depth, density, temperature, &
    & age, grain_radius) result(d_density)
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

    integer, intent(in)                   :: nz
    real(prec), intent(in)                :: dt
    real(prec), dimension(nz), intent(in) :: depth
    real(prec), dimension(nz), intent(in) :: density
    real(prec), dimension(nz), intent(in) :: temperature
    real(prec), dimension(nz), intent(in) :: age
    real(prec), dimension(nz), intent(in) :: grain_radius

    real(prec), dimension(nz) :: d_density
    real(prec)                :: mean_temperature
    integer                   :: n

    call tfm_density_do_nothing(nz, grain_radius)

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
    &  (/ a0, ALPHA0, EC0 /), &
    &  (/ a1, ALPHA1, EC1 /), &
    &  depth,                 &
    &  temperature,           &
    &  density,               &
    &  age,                   &
    &  d_density              &
    )
  end function tfm_density_arthern2010


  function tfm_density_ligtenberg2011(nz, dt, depth, density, temperature, &
    & age, grain_radius) result(d_density)
    implicit none

    integer, intent(in)                   :: nz
    real(prec), intent(in)                :: dt
    real(prec), dimension(nz), intent(in) :: depth
    real(prec), dimension(nz), intent(in) :: density
    real(prec), dimension(nz), intent(in) :: temperature
    real(prec), dimension(nz), intent(in) :: age
    real(prec), dimension(nz), intent(in) :: grain_radius

    real(prec), dimension(nz) :: d_density
    real(prec), dimension(nz) :: mean_acc
    integer                   :: n

    d_density = tfm_density_arthern2010(nz, dt, depth, density, &
      & temperature, age, grain_radius)

    ! boundary between first and seconds stage
    do n = nz, 1, -1
      if ( density(n) >= 550.0 ) EXIT
    end do

    ! change accumulation to (kg a-1 m-2)
    call tfm_density_mean_acc(nz, depth, density, age, mean_acc)
    mean_acc = mean_acc * WATER_DENSITY

    d_density(n+1:nz-1) = (                             &
    &  d_density(n+1:nz-1)                              &
    &  * (1.435 - (0.151 * log(mean_acc(n+1:nz-1))))    &
    )
    d_density(1:n) = (                              &
    &  d_density(1:n)                               &
    &  * (2.366 - (0.293 * log(mean_acc(1:n))))     &
    )
  end function tfm_density_ligtenberg2011


  function tfm_density_simonsen2013(nz, dt, depth, density, temperature, &
    & age, grain_radius) result(d_density)
    implicit none

    ! fist stage parameters (as implemented in CFM)
    real(prec), parameter :: F0 = 0.8

    ! second stage parameters (as implemented in CFm)
    real(prec), parameter :: F1 = 1.25

    integer, intent(in)                   :: nz
    real(prec), intent(in)                :: dt
    real(prec), dimension(nz), intent(in) :: depth
    real(prec), dimension(nz), intent(in) :: density
    real(prec), dimension(nz), intent(in) :: temperature
    real(prec), dimension(nz), intent(in) :: age
    real(prec), dimension(nz), intent(in) :: grain_radius

    real(prec), dimension(nz) :: d_density

    integer                   :: n
    real(prec), dimension(nz) :: mean_acc
    real(prec)                :: mean_temperature

    d_density = tfm_density_arthern2010(nz, dt, depth, density, &
      & temperature, age, grain_radius)

    ! boundary between first and seconds stage
    do n = nz, 1, -1
      if ( density(n) >= 550.0 ) EXIT
    end do

    ! change accumulation to (kg a-1 m-2)
    call tfm_density_mean_acc(nz, depth, density, age, mean_acc)
    mean_acc = mean_acc * WATER_DENSITY

    ! 10 m temperature
    do n = nz, 1, -1
      if ( (depth(n) - depth(nz)) <= -10.0 ) EXIT
    end do
    mean_temperature = temperature(n)

    d_density(n+1:nz-1) = F0 * d_density(n+1:nz-1)
    d_density(1:n) = (                                 &
    &  d_density(1:n)                                  &
    &  * (61.7 / (mean_acc(1:n)**0.5))                 &
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
