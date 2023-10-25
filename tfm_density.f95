module tfm_density_tools
  use tfm_essentials
  use tfm_constants
  implicit none
! ---------------------------------------------------------------------
! Module: tfm_density_tools
!
! This module contains procedures commonly used by density related
! functions and subroutines.
!
! Dependencies: tfm_essentials, tfm_constants
!
! Interfaces:
!   tfm_density_arrhenius: Interface for functions
!                          "tfm_density_arrhenius_i" and
!                          "tfm_density_arrhenius_n".
!
! Functions:
!   tfm_density_arrhenius_i: Arrhenius equation for dimension 1.
!   tfm_density_arrhenius_n: Arrhenius equation for dimension n.
!
! Subroutines:
!   tfm_density_lin_interp: Simple linear interpolation.
!   tfm_density_bagmean: Bagmean from firn profile.
! ---------------------------------------------------------------------

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
! ---------------------------------------------------------------------
! Function: tfm_density_arrhenius_i
!
! The function computes an arrhenius equation for a given factor, 
! activation energy and temperature. The dimension of the input and
! output is always one.
!
! Author: Timm Schultz
! 
! Arguments:
!   i: Dummy variable needed for the interface (tfm_density_arrhenius).
!   factor: Pre factor of the arrhenius equation.
!   activation_energy: Activation energy of the arrhenius equation.
!   tempreature: Temperature (K).
!
! Result:
!   arrhenius: Arrhenius Factor (dimension 1).
! ---------------------------------------------------------------------

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
! ---------------------------------------------------------------------
! Function: tfm_density_arrhenius_n
!
! The function computes an arrhenius equation for a given factor, 
! activation energy and temperature. The dimension of the input and
! output is defined by input.
!
! Author: Timm Schultz
! 
! Arguments:
!   n: Dimension of temperature input and arrhenius factor output.
!   factor: Pre factor of the arrhenius equation.
!   activation_energy: Activation energy of the arrhenius equation.
!   temperature: Temperature (K).
!
! Result:
!   arrhenius: Arrhenius factor (dimension n).
! ---------------------------------------------------------------------

    arrhenius = factor * exp(-activation_energy / (GAS_CONST * temperature))
  end function tfm_density_arrhenius_n


  subroutine tfm_density_lin_interp(z0, z1, v0, v1, dz, v)
    implicit none

    real(prec), intent(in) :: z0, z1
    real(prec), intent(in) :: v0, v1
    real(prec), intent(in) :: dz

    real(prec), intent(inout) :: v
! ---------------------------------------------------------------------
! Subroutine: tfm_density_lin_interp
!
! Simple routine for linear interpolation between given values.
!
! Author: Timm Schultz
!
! Arguments:
!   z0: First given x-value.
!   z1: Second given x-value.
!   v0: First fiven y-value.
!   v1: Second given y-value.
!   dz: Distance to value z0 at which to interpolate at.
!   v - on input: Variable to store the result.
!
! Result:
!   v - on output: Value from linear interpolation at z0 + dz.
! ---------------------------------------------------------------------

    v = v0 + ((v1 - v0) / (z1 - z0)) * dz
  end subroutine tfm_density_lin_interp


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
! ---------------------------------------------------------------------
! Subroutine: tfm_density_bagmean
!
! The routine computes bagmean values of a given density profile.
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of arrays "depth" and "density".
!   depth: Array containing the depth values of the density profile (m).
!   density: Array containing the density values of the density
!            profile (kg m**-3).
!   dz: Size of the "window" the bagmean values are computed from (m).
!   mz: Length of the bagmean density profile (total length of
!       the profile is mz*dz)
!   bagmean - on input: Variable to store the bagmean profile.
!
! Result:
!   bagmean - on output: Array of dimension "mz" contaning the depth
!                        (1,:) values (m) of the bagmean profile and
!                        the corresponding density (2,:)
!                        values (kg m**-3).
! ---------------------------------------------------------------------

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
  use tfm_essentials
  use tfm_constants
  implicit none
! ---------------------------------------------------------------------
! Module: tfm_density_stress
!
! The module contains common stress related procedures.
!
! Dependencies: tfm_essentials, tfm_constants
!
! Functions:
!  tfm_density_boyleMariotte: Boyle-Mariotte factor from absolute
!                             density.
!  tfm_rel_density_boyleMariotte: Boyle-Mariotte factor from
!                                 relative density
!
! Subroutines:
!  tfm_density_computeStress: stres from depth and density
! ---------------------------------------------------------------------

  contains


  subroutine tfm_density_computeStress(nz, depth, density, stress)
    implicit none

    integer, intent(in)                      :: nz
    real(prec), dimension(nz), intent(in)    :: depth
    real(prec), dimension(nz), intent(in)    :: density
    real(prec), dimension(nz), intent(inout) :: stress

    integer                   :: n
    real(prec), dimension(nz) :: dz
! ---------------------------------------------------------------------
! Subroutine tfm_density_computeStress
!
! Routine to compute the stress along a given firn profile, by
! integration of the density.
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of variables "depth", "density", and "stress".
!   depth: Depth along the firn profile (m).
!   density: Density along the firn profile (kg m**-3).
!   stress - on input: Variable to store the result.
!
! Result:
!   stress - on output: Stress computed from "depth" and
!                       "density" (Pa).
! ---------------------------------------------------------------------

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
! ---------------------------------------------------------------------
! Function: tfm_density_boyleMariotte
!
! Boyle-Mariotte equation to compute the pore pressure within bubbly
! ice from given absoulte density.
!
! Author: Timm Schultz
!
! Arguments:
!   density: Absolute density (kg m**-3).
!
! Result:
!   boyle_mariotte: Boyle-Mariotte factor. The pore pressure is the
!                   product of this factor and the pressure at pore
!                   close off (1).
! ---------------------------------------------------------------------

    boyle_mariotte = (                                &
    &  (density * (ICE_DENSITY - CLOSEOFF_DENSITY))   &
    &  / (CLOSEOFF_DENSITY * (ICE_DENSITY - density)) &
    )
  end function tfm_density_boyleMariotte


  function tfm_density_rel_boyleMariotte(rel_density) result(rel_boyle_mariotte)
    implicit none

    real(prec), intent(in) :: rel_density
    real(prec)             :: rel_boyle_mariotte
! ---------------------------------------------------------------------
! Function: tfm_density_rel_boyleMariotte
!
! Boyle-Mariotte equation to compute the pore pressure within bubbly
! ice from given relative density.
!
! Author: Timm Schultz
!
! Arguments:
!   density: Relative density (1).
!
! Result:
!   boyle_mariotte: Boyle-Mariotte factor. The pore pressure is the
!   product of this factor and the pressure at pore close off (1).
! ---------------------------------------------------------------------

    rel_boyle_mariotte = (                                             &
    &  rel_density * (1.0_prec - (CLOSEOFF_DENSITY / ICE_DENSITY))     &
    &  / ((CLOSEOFF_DENSITY / ICE_DENSITY) * (1.0_prec - rel_density)) &
    )
  end function tfm_density_rel_boyleMariotte
end module tfm_density_stress



module tfm_density_processes
  use tfm_essentials
  use tfm_constants
  use tfm_density_tools
  implicit none

! ---------------------------------------------------------------------
! Module: tfm_density_processes
!
! The module contains functions for process based firn densification
! modelling.
!
! Dependencies: tfm_essentials, tfm_constants, tfm_density_tools
!
! Functions:
!   tfm_density_NabarroHerringCreep
!   tfm_density_NabarroHerringCreepMod
!   tfm_density_CobleCreep
!   tfm_density_DislocationCreep
!   tfm_density_DislocationCreepMod
!   tfm_density_GrainBoundarySliding
! ---------------------------------------------------------------------

  contains

  function tfm_density_NabarroHerringCreep(temperature, density, grain_radius, &
    & driving_force) result(strain_rate)
    implicit none

    real(prec), intent(in) :: temperature
    real(prec), intent(in) :: density
    real(prec), intent(in) :: grain_radius
    real(prec), intent(in) :: driving_force
    real(prec)             :: strain_rate

    real(prec) :: rate_factor_lattice_diffusion
! ---------------------------------------------------------------------
! Function: tfm_density_NabarroHerringCreep
!
! Maeno, N. and Ebinuma, T. Pressure Sintering of Ice and Its
! Implication to the Densification of Snow at Polar Glaciers and Ice
! Sheets. J. Phys. Chem., 87, 4103-4110, (1983). 
!
! Author: Timm Schultz
!
! Arguments:
!   temperature: Temperaure (K).
!   grain_radius: Grain radius (m).
!   driving force: Sintering driving force (Pa).
!
! Result:
!   strain_rate: Resulting strain rate (s**-1).
! ---------------------------------------------------------------------

    rate_factor_lattice_diffusion = tfm_density_arrhenius( &
    &  1,                                                  &
    &  PRE_FACTOR_LATTICE_DIFFUSION,                       &
    &  ACTIVATION_ENERGY_LATTICE_DIFFUSION,                &
    &  temperature                                         &
    )
    
    !strain_rate = (                                                       &
    !&  -((10.0_prec * VOLUME_H2O) / (3.0_prec * BOLTZMANN * temperature)) &
    !&  * (1.0_prec / (grain_radius**2.0_prec))                            &
    !&  * rate_factor_lattice_diffusion                                    &
    !&  * driving_force                                                    &
    !)

    strain_rate = (                                                       &
    &  -((10.0_prec * VOLUME_H2O) / (2.0_prec * BOLTZMANN * temperature)) &
    &  * (1.0_prec / (grain_radius**3.0_prec))                            &
    &  * rate_factor_lattice_diffusion                                    &
    &  * (ICE_DENSITY / density) &
    &  * driving_force                                                    &
    )
  end function tfm_density_NabarroHerringCreep


  function tfm_density_NabarroHerringCreepMod(temperature, density, &
    & grain_radius, pore_radius, driving_force) result(strain_rate)
    implicit none

    real(prec), intent(in) :: temperature
    real(prec), intent(in) :: density
    real(prec), intent(in) :: grain_radius
    real(prec), intent(in) :: pore_radius
    real(prec), intent(in) :: driving_force
    real(prec)             :: strain_rate

    real(prec) :: rate_factor_lattice_diffusion
! ---------------------------------------------------------------------
! Function: tfm_density_NabarroHerringCreepMod
!
! Maeno, N. and Ebinuma, T. Pressure Sintering of Ice and Its
! Implication to the Densification of Snow at Polar Glaciers and Ice
! Sheets. J. Phys. Chem., 87, 4103-4110, (1983). 
!
! Author: Timm Schultz
!
! Arguments:
!   temperature: Temperature (K).
!   density: Density (kg m**-3).
!   grain_radius: Grain radius (m).
!   pore_radius: Pore radius (m).
!   driving_force: Sintering dirving force (Pa).
!
! Result:
!   strain_rate: Resulting strain rate (s**-1).
! ---------------------------------------------------------------------
    
    rate_factor_lattice_diffusion = tfm_density_arrhenius( &
    &  1,                                                  &
    &  PRE_FACTOR_LATTICE_DIFFUSION,                       &
    &  ACTIVATION_ENERGY_LATTICE_DIFFUSION,                &
    &  temperature                                         &
    )

    strain_rate = (                                                  &
    &  -((3.0_prec * VOLUME_H2O) / (BOLTZMANN * temperature))        &
    &  * (1.0_prec / (grain_radius**2.0_prec))                       &
    &  * (density / ICE_DENSITY)                                     &
    &  * ((1.0_prec / grain_radius) * (                              &
    &    (pore_radius * grain_radius) / (grain_radius - pore_radius) &
    &  ))                                                            &
    &  * rate_factor_lattice_diffusion                               &
    &  * driving_force                                               &
    )
  end function tfm_density_NabarroHerringCreepMod


  function tfm_density_CobleCreep(temperature, grain_radius, &
    & driving_force) result(strain_rate)
    implicit none
    
    real(prec), intent(in) :: temperature
    real(prec), intent(in) :: grain_radius
    real(prec), intent(in) :: driving_force
    real(prec)             :: strain_rate

    real(prec) :: rate_factor_boundary_diffusion
! ---------------------------------------------------------------------
! Function: tfm_density_CobleCreep
!
! Maeno, N. and Ebinuma, T. Pressure Sintering of Ice and Its
! Implication to the Densification of Snow at Polar Glaciers and Ice
! Sheets. J. Phys. Chem., 87, 4103-4110, (1983). 
!
! Author: Timm Schultz
!
! Arguments:
!   temperature: Temeprature (K).
!   grain_radius: Grain raidus (m).
!   driving_force: Sintering driving force (Pa).
!
! Result:
!   strain_rate: Resulting strain rate (s**-1).
! ---------------------------------------------------------------------
    
    rate_factor_boundary_diffusion = tfm_density_arrhenius( &
    &  1,                                                   &
    &  PRE_FACTOR_BOUNDARY_DIFFUSION,                       &
    &  ACTIVATION_ENERGY_BOUNDARY_DIFFUSION,                &
    &  temperature                                          &
    )

    strain_rate = (                                           &
    &  -(                                                     &
    &    (37.0_prec * VOLUME_H2O * 2.0_prec * BURGERS_VECTOR) &
    &    / (2.0_prec * BOLTZMANN * temperature)               &
    &  )                                                      &
    &  * (1.0_prec / (grain_radius**3.0_prec))                &
    &  * rate_factor_boundary_diffusion                       &
    &  * driving_force                                        &
    )
  end function tfm_density_CobleCreep


  function tfm_density_DislocationCreep(temperature, density, &
    & driving_force) result(strain_rate)
    implicit none

    real(prec), intent(in) :: temperature
    real(prec), intent(in) :: density
    real(prec), intent(in) :: driving_force
    real(prec)             :: strain_rate

    real(prec) :: rate_factor_dislocation_creep
! ---------------------------------------------------------------------
! Function: tfm_density_DislocationCreep
!
! Maeno, N. and Ebinuma, T. Pressure Sintering of Ice and Its
! Implication to the Densification of Snow at Polar Glaciers and Ice
! Sheets. J. Phys. Chem., 87, 4103-4110, (1983). 
!
! Author: Timm Schultz
!
! Arguments:
!   temperature: Temeprature (K).
!   density: Density (kg m**-3).
!   driving_force: Sintering driving force (Pa).
!
! Result:
!   strain_rate: Resulting strain rate (s**-1).
! ---------------------------------------------------------------------
    
    if (temperature > 263.15_prec) then
      rate_factor_dislocation_creep = tfm_density_arrhenius( &
      &  1,                                                  &
      &  PRE_FACTOR_DISLOCATION_CREEP_HIGH,                  &
      &  ACTIVATION_ENERGY_DISLOCATION_CREEP_HIGH,           &
      &  temperature                                         &
      )

    else if (temperature <= 263.15) then
      rate_factor_dislocation_creep = tfm_density_arrhenius( &
      &  1,                                                  &
      &  PRE_FACTOR_DISLOCATION_CREEP_LOW,                   &
      &  ACTIVATION_ENERGY_DISLOCATION_CREEP_LOW,            &
      &  temperature                                         &
      )

    else
       print *, 'Module: tfm_density_processes'
       print *, 'Function: tfm_density_DislocationCreep'
       print *, ''
       print *, 'The temperature seems to show and irregular value.'
       print *, 'Stopping right here!'
       STOP
    end if

    strain_rate = (                                                                         &
    &  (-2.0_prec * rate_factor_dislocation_creep)                                          &
    &  * (                                                                                  &
    &    (1.0_prec - (density / ICE_DENSITY))                                               &
    &    / ((1.0_prec - ((1.0_prec - (density / ICE_DENSITY))**(1.0_prec / ICE_N)))**ICE_N) &
    &  )                                                                                    &
    &  * (((2.0_prec / ICE_N) * driving_force)**ICE_N)                                      &
    )
  end function tfm_density_DislocationCreep


  function tfm_density_DislocationCreepMod(temperature, density, &
    & driving_force) result(strain_rate)
    implicit none

    real(prec), intent(in) :: temperature
    real(prec), intent(in) :: density
    real(prec), intent(in) :: driving_force
    real(prec)             :: strain_rate

    real(prec) :: rate_factor_dislocation_creep
! ---------------------------------------------------------------------
! Function: tfm_density_DislocationCreep
!
! Maeno, N. and Ebinuma, T. Pressure Sintering of Ice and Its
! Implication to the Densification of Snow at Polar Glaciers and Ice
! Sheets. J. Phys. Chem., 87, 4103-4110, (1983). 
!
! Author: Timm Schultz
!
! Arguments:
!   temperature: Temeprature (K).
!   density: Density (kg m**-3).
!   driving_force: Sintering driving force (Pa).
!
! Result:
!   strain_rate: Resulting strain rate (s**-1).
! ---------------------------------------------------------------------
    
    if (temperature > 263.15_prec) then
      rate_factor_dislocation_creep = tfm_density_arrhenius( &
      &  1,                                                  &
      &  PRE_FACTOR_DISLOCATION_CREEP_HIGH,                  &
      &  ACTIVATION_ENERGY_DISLOCATION_CREEP_HIGH,           &
      &  temperature                                         &
      )

    else if (temperature <= 263.15) then
      rate_factor_dislocation_creep = tfm_density_arrhenius( &
      &  1,                                                  &
      &  PRE_FACTOR_DISLOCATION_CREEP_LOW,                   &
      &  ACTIVATION_ENERGY_DISLOCATION_CREEP_LOW,            &
      &  temperature                                         &
      )

    else
       print *, 'Module: tfm_density_processes'
       print *, 'Function: tfm_density_DislocationCreepMod'
       print *, ''
       print *, 'The temperature seems to show and irregular value.'
       print *, 'Stopping right here!'
       STOP
    end if

    strain_rate = (                                                                         &
    &  (-3.0_prec / 2.0_prec)                                                               &
    &  * rate_factor_dislocation_creep                                                      &
    &  * (                                                                                  &
    &    (1.0_prec - (density / ICE_DENSITY))                                               &
    &    / ((1.0_prec - ((1.0_prec - (density / ICE_DENSITY))**(1.0_prec / ICE_N)))**ICE_N) &
    &  )                                                                                    &
    &  * (((3.0_prec / (2.0_prec * ICE_N)) * driving_force)**ICE_N)                         &
    )
  end function tfm_density_DislocationCreepMod


  function tfm_density_GrainBoundarySliding(temperature, density, &
    & grain_radius, neck_radius, amplitude_gbo, driving_force)    &
    & result(strain_rate)
    implicit none
    
    real(prec), intent(in) :: temperature
    real(prec), intent(in) :: density
    real(prec), intent(in) :: grain_radius
    real(prec), intent(in) :: neck_radius
    real(prec), intent(in) :: amplitude_gbo
    real(prec), intent(in) :: driving_force
    real(prec)             :: strain_rate

    real(prec) :: rate_factor_boundary_diffusion
! ---------------------------------------------------------------------
! Function: tfm_density_GrainBoundarySliding
!
! Alley, R. B. Firn Densification by Grain Boundary Sliding: A First
! Model. J. Phys. Colloques, 48, C1, C1-249-C1-256, (1987), https://
! doi.org/10.1051/jphyscol:1987135
!
! Author: Timm Schultz
!
! Arguments:
!   temperature: Temperature (K).
!   density: Density (kg m**-3).
!   grain_radius: Grain radius (m).
!   neck_radius: Radius of necks between grains (m).
!   amplitude_gbo: Amplitude of grain boundary obstructions (m).
!   driving_force: Sintering driving force (Pa).
!
! Reuslt:
!   strain_rate: Resulting strain rate (s**-1).
! ---------------------------------------------------------------------
    
    rate_factor_boundary_diffusion = tfm_density_arrhenius( &
    &  1,                                                   &
    &  PRE_FACTOR_BOUNDARY_DIFFUSION,                       &
    &  ACTIVATION_ENERGY_BOUNDARY_DIFFUSION,                &
    &  temperature                                          &
    )

    strain_rate = (                                                     &
    &  (-2.0_prec / 15.0_prec)                                          &
    &  * (2.0_prec * BURGERS_VECTOR)                                    &
    &  * (                                                              &
    &    (8.0_prec * rate_factor_boundary_diffusion * VOLUME_H2O)       &
    &    / (BOLTZMANN * temperature * (amplitude_gbo**2.0_prec))        &
    &  )                                                                &
    &  * (grain_radius / (neck_radius**2.0_prec))                       &
    &  * ((ICE_DENSITY / density)**3.0_prec)                            &
    &  * (1.0_prec - ((5.0_prec * density) / (3.0_prec * ICE_DENSITY))) &
    &  * driving_force                                                  &
    )
  end function tfm_density_GrainBoundarySliding
end module tfm_density_processes



module tfm_density_herronLangway
  use tfm_essentials 
  use tfm_constants
  implicit none
! ---------------------------------------------------------------------
! Module: tfm_density_herronLangway
!
! This module contains generic routines to solve firn densification
! models of the Herron & Langway type. These routines can be used to
! solve various firn densification models based on this model type.
!
! Herron, M. M. and Langway, C. C. Firn Densification: An Empirical
! Model. Journal of Glaciology, 25 (93), 373-385, (1980). https://
! doi.org/10.3189/S0022143000015239
!
! Dependencies: tfm_essentials, tfm_constants
!
! Subroutines:
!   tfm_density_HLtype: Generic function for Herron & Langway type
!                       models.
! ---------------------------------------------------------------------

  contains


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
! ---------------------------------------------------------------------
! Subroutine: tfm_density_HLtype
!
! The routine computes the density change resulting from a "generic"
! densification model of the Herron & Langway type.
!
! Herron, M. M. and Langway, C. C. Firn Densification: An Empirical
! Model. Journal of Glaciology, 25 (93), 373-385, (1980). https://
! doi.org/10.3189/S0022143000015239
!
! Many firn densification models rely on this type of model. The idea
! is to formulate a generic function following the general form of the
! material model. It is then solved by passing three arguments for the
! first and the second stage of firn densification. All models following
! the concept of Herron & Langway can be broken down to these six
! parameters.
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension variables "depth", "temperature", "density", "age",
!       "d_density".
!   dt: Time step (s).
!   stage1_params: Parameters for the first densification stage
!                  (dimension (3)).
!   stage2_params: Parameters for the seconds densification stage
!                  (dimension (3)).
!   depth: Depth of the firn profile (m).
!   temperature: Temperature of the firn profile (K).
!   density: Density of the firn profile (kg m**-3).
!   age: Age of the firn profile (s).
!   d_density - on input: Variable for the result.
!
! Result:
!   d_density - on output: Change of the density (kg m**-3).
! ---------------------------------------------------------------------

    ! computation of the mean accumulation rate 
    ! over the lifetime of the firn parcel (kg a-1 m-2)
    call tfm_essentials_mean_acc(nz, depth, density, age, mean_acc)
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
  use tfm_essentials
  use tfm_constants
  implicit none

  interface tfm_density_arztCoordination
    module procedure tfm_density_arztCoordination_i, tfm_density_arztCoordination_n
  end interface tfm_density_arztCoordination

  interface tfm_density_arztContactarea
    module procedure tfm_density_arztContactarea_i, tfm_density_arztContactarea_n
  end interface tfm_density_arztContactarea

  ! parameters
  real(prec), parameter :: ZZERO = 7.3_prec
  !real(prec), parameter :: ZZERO = 4.81_prec
  real(prec), parameter :: CARZT = 15.5_prec
! ---------------------------------------------------------------------
! Module: tfm_density_fischmeisterArzt
!
! The module contains routines to compute the coordination number and
! contact area of particles during sintering based on the density
! following the work of Fischmeister & Arzt. For different uses
! different functions handling different dimensions are connected via
! suitable interfaces.
!
! Fischmeister, H. F. and Arzt, E. Densification of powders by particle
! deformation. Powder Metallurgy, 26 (2), 82-88, (1983).
! https://doi.org/10.1179/pom.1983.26.2.82
!
! Dependencies: tfm_essentials, tfm_constants
! 
! Interfaces:
!   tfm_density_arztCoordination: Interface for coordination number
!                                 ("tfm_density_arztCoordination_i",
!                                 "tfm_arztCoordination_n").
!   tfm_density_arztContactarea: Interface for contact area
!                                ("tfm_density_arztContactarea_i",
!                                "tfm_arztContactarea_n").
!
! Parameters:
!   ZZERO: Coordination number at densest packing (1).
!   CARZT: Factor (1)
!
! Functions:
!   tfm_density_arztCoordination_i: Coordination number, dimension (1).
!   tfm_density_arztCoordination_n: Coordination number, dimenison (n).
!   tfm_density_arztContactarea_i: Contact area, dimension(1).
!   tfm_density_arztContactarea_n: Contact area, dimension(n).
! ---------------------------------------------------------------------

  contains


  function tfm_density_arztCoordination_i(n, rel_density, d_zero) &
    & result(coordination_number)
    implicit none

    integer, intent(in)    :: n
    real(prec), intent(in) :: rel_density
    real(prec), intent(in) :: d_zero
    real(prec)             :: coordination_number
    integer                :: m
! ---------------------------------------------------------------------
! Function tfm_density_arztCoordination_i
!
! Computation of the coordination number following Fischmeister & Arzt
! (1983) for dimension (1).
!
! Author: Timm Schultz
!
! Arguments:
!   n: Dummy argument for the interface (always = 1).
!   rel_density: Relative density (1).
!   d_zero: Relative density at densest packing (1).
!
! Result:
!   coordination_number: Coordination number at given density (1).
! ---------------------------------------------------------------------

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
! ---------------------------------------------------------------------
! Function tfm_density_arztCoordination_n
!
! Computation of the coordination number following Fischmeister & Arzt
! (1983) for dimension (n).
!
! Author: Timm Schultz
!
! Arguments:
!   n: Dimension of variables "rel_density", "coordination_number".
!   rel_density: Relative density (1).
!   d_zero: Relative density at densest packing (1).
!
! Result:
!   coordination_number: Coordination number at given density (1).
! ---------------------------------------------------------------------

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
! ---------------------------------------------------------------------
! Function tfm_density_arztContactarea_i
!
! Computation of the contact area between grains  following
! Fischmeister & Arzt (1983) for dimension (1).
!
! Author: Timm Schultz
!
! Arguments:
!   n: Dummy argument for the interface (always = 1).
!   rel_density: Relative density (1).
!   d_zero: Relative density at densest packing (1).
!
! Result:
!   contact_area: Contact area at given density (1).
! ---------------------------------------------------------------------

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
! ---------------------------------------------------------------------
! Function tfm_density_arztContactarea_n
!
! Computation of the contact area between grains  following
! Fischmeister & Arzt (1983) for dimension (n).
!
! Author: Timm Schultz
!
! Arguments:
!   n: Dimension of variables "rel_density", "coordination_number".
!   rel_density: Relative density (1).
!   d_zero: Relative density at densest packing (1).
!
! Result:
!   contact_area: Contact area at given density (1).
! ---------------------------------------------------------------------

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
  use tfm_essentials
  use tfm_constants
  use tfm_density_tools

  real(prec), parameter :: STAGE_DIV = 0.81_prec

  interface
    function invariant_inter(nz, param_a, param_b, strain_rate_inp) &
      & result(invariant)
      use tfm_essentials
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
      use tfm_essentials
      implicit none

      integer, intent(in)                   :: nz
      real(prec), dimension(nz), intent(in) :: param
      real(prec), dimension(nz), intent(in) :: rate_factor
      real(prec), dimension(nz), intent(in) :: invariant

      real(prec), dimension(nz) :: viscosity
    end function viscosity_inter
  end interface
! ---------------------------------------------------------------------
! Module: tfm_density_gagliardini
!
! The module contains routines to solve a firn densification model of
! the kind first described for firn by:
!
! Gagliardini, O. and Meyssonnier, J. Flow simulation of a firn-covered
! cold glacier. Annals of Glaciology, 24, 242-248 (1997).
! https://doi.org/10.3189/S0260305500012246
!
! The model is based on the work of Duva & Crow:
!
! Duva, J. M. and Crow, P. D. Analysis of consolidation of reinforced
! materials by power-law creep. Mechanics of Materials, 17 (1), 25-32,
! (1994). https://doi.org/10.1016/0167-6636(94)90011-6
!
! Dependencies: tfm_essentials, tfm_constants, tfm_denisty_tools
! 
! Interfaces:
!   invariant_inter: Interface for the invariant function.
!   viscosity_inter: Interface for viscosity functions.
!
! Parameters:
!   STAGE_DIV: Relative density at transition stage II - stage III (1).
! 
! Functions:
!   tfm_density_gagliardiniParamA0: Parameter A0 of the model.
!   tfm_density_gagliardiniParamB0: Parameter B0 of the model.
!   tfm_density_gagliardiniRate: Rate factor.
!   tfm_density_gagliradiniSolve: Function for solving the model
!                                 interatively.
! ---------------------------------------------------------------------

  contains


  function tfm_density_gagliardiniParamA0(density) result(param_a0)
    implicit none

    real(prec), intent(in) :: density
    real(prec)             :: param_a0
    real(prec)             :: rel_density
! ---------------------------------------------------------------------
! Function: tfm_density_gagliardiniParmA0
!
! Parameter a0 of the firn densification model as described by
! Gagliardini & Meyssonnier (1997), following Duva & Crow (1994).
!
! Author: Timm Schultz
!
! Arguments:
!   density: Absolute density (kg m**-3).
!
! Result:
!   param_a0: Parameter a0.
! ---------------------------------------------------------------------

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
! ---------------------------------------------------------------------
! Function: tfm_density_gagliardiniParamB0
!
! Parameter b0 of the firn densification model as described by
! Gagliardini & Meyssonnier (1997), following Duva & Crow (1994).
!
! Author: Timm Schultz
!
! Arguments:
!   density: Absolute density (kg m**-3).
!
! Result:
!   param_b0: Parameter b0.
! ---------------------------------------------------------------------
    
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
! ---------------------------------------------------------------------
! Function: tfm_density_gagliardiniRate
!
! Rate factor for ice flow following:
! Greve, R. and Blatter, H. Dynamics of Ice Sheets and Galciers.
! Springer, Berlin, (2009).
!
! Author: Timm Schultz
!
! Parameters (Greve & Blatter, 2009):
!   TEMP_DIV: Temperature dividing high and low temperature (K).
!   PRE_FACTOR_LOW: Pre factor of the arrhenius equation at low
!                   temperatures.
!   PRE_FACTOR_HIGH: Pre factor of the arrhenius equation at high
!                    temperatures.
!   ACTIVATION_ENERGY_LOW: Activation energy of the arrhenius equation
!                          at low temperatures.
!   ACTIVATION_ENERGY_HIGH: Activation energy of the arrhenius equation
!                           at hight temperatures.
!
! Arguments:
!   nz: Dimension of variables "temperature" and "rate_factor".
!   temperature: Temperature (K).
!
! Result:
!   rate_factor: Rate factor.
! ---------------------------------------------------------------------
    
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
! ---------------------------------------------------------------------
! Function tfm_density_gagliardiniSolve
!
! This function solves a firn densification model of the kind described
! by Gagliardini & Meyssonnier (1997). Because this kind of model is
! non-linear it is solved iteratively.
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of variables "density", "stress", "param_a",
!       "param_b", "rate_factor".
!   density: Absolute density (kg m**-3).
!   stress: Stress due to overburden firn (Pa).
!   dt: Time step (s).
!   param_a: Parameter a along the firn profile.
!   param_b: Parameter b along the firn profile.
!   rate_factor: Ice flow rate factor.
!   invariant_func: Function defining the strain invariant (procedure,
!                   interface: invariant_func).
!   shear_visco_func: Function defining the shear viscosity (procedure,
!                     interface: viscosity_inter).
!   bulk_visco_func: Function defining the bulk viscosity (procedure,
!                    interface: viscosity_inter).
!
! Result:
!   d_density: Density change along the firn profile (kg m**-3).
! ---------------------------------------------------------------------

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
      invariant = invariant_func(          &
      &  nz, param_a, param_b, strain_rate &
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
      if ( (density(n) + d_density(n)) > (ICE_DENSITY - 10.0e-5_prec) ) then
        d_density(n) = (ICE_DENSITY - density(n)) - 10.0e-5_prec
      end if
    end do
  end function tfm_density_gagliardiniSolve
end module tfm_density_gagliardini



module tfm_density
  use tfm_essentials
  use tfm_constants
  use tfm_density_tools
  use tfm_density_herronLangway
  use tfm_density_stress
  use tfm_density_fischmeisterArzt
  use tfm_density_processes
  use tfm_density_gagliardini
  implicit none
! ---------------------------------------------------------------------
! Module: tfm_density
!
! This module is a collection of various different firn densification
! models. All functions compute the change in density. The arguments
! passed to the individual functions are always the same allowing to
! pass them via a suitable interface. This means that some arguments
! may not be used by the function.
!
! Dependencies: tfm_essentials, tfm_constants, tfm_density_tools,
!               tfm_density_herronLangway, tfm_density_stress,
!               tfm_density_fischmeisterArzt, tfm_density_processes,
!               tfm_density_gagliardini
!
! Functions:
!   tfm_density_depth: Change in depth due to densification.
!   tfm_density_gagliardini1998: Gagliardini & Meyssonier (1998).
!   tfm_density_timmsfit: Timm's version of the Gagliardini model.
!   tfm_density_greve2009: Greve & Blatter (2009).
!   tfm_density_zwinger2007: Zwinger et al. (2007).
!   tfm_density_breant2017: Breant et al. (2017).
!   tfm_density_medley2020: Medley et al. (2020) (preprint).
!   tfm_density_herron1980: Herron & Langway (1980) (transient).
!   tfm_density_arthern1998: Arthern & Wingham (1998).
!   tfm_density_li2003: Li & Zwally (2003).
!   tfm_density_helsen2008: Helsen et al. (2008).
!   tfm_density_arthern2010: Arthern et al. (2010).
!   tfm_density_ligtenberg2011: Ligtenberg et al. (2011).
!   tfm_density_simonsen2013: Simonsen et al. (2013).
!
! Subroutines:
!   tfm_density_herron1980_analytical: Herron & Langway (analytical).
! ---------------------------------------------------------------------

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
! ---------------------------------------------------------------------
! Function: tfm_density_depth
!
! The function computes the change in depth along a firn profile from
! change in density.
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of variables "depth", "density", "d_density",
!       "d_depth".
!   depth: Depth along the firn profile (m).
!   density: Density along the firn profile (kg m**-3).
!   d_density: Density change along the firn profile (kg m**-3).
!
! Result:
!   d_depth: Change in depth along the firn profile (m).
! ---------------------------------------------------------------------

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
! ----------------------------------------------------------------------
! Function: tfm_denisty_gagliardini1998
!
! Gagliardini, O. and Meyssonnier, J. Flow simulation of a firn-covered
! cold glacier. Annals of Glaciology, 24, 242-248, (1997).
! https://doi.org/10.3189/S0260305500012246
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of variables "depth", "density", "temperature", "age",
!       "grain_radius", "d_density".
!   dt: Time step (s).
!   depth: Depth along the firn profile (m).
!   density: Density along the firn profile (kg m**-3).
!   temperature: Temeperature along the firn profile (K).
!   age: Age along the firn profile (s).
!   grain_radius: Grain radius along the firn profile (m).
!
! Result:
!   d_density: Density change along the firn profile (kg m**-3).
! ----------------------------------------------------------------------

    ! doing nothing
    call tfm_essentials_do_nothing(nz, age)
    call tfm_essentials_do_nothing(nz, grain_radius)

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
! ----------------------------------------------------------------------
! Function: tfm_density_timmsfit
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of variables "depth", "density", "temperature", "age",
!       "grain_radius", "d_density".
!   dt: Time step (s).
!   depth: Depth along the firn profile (m).
!   density: Density along the firn profile (kg m**-3).
!   temperature: Temeperature along the firn profile (K).
!   age: Age along the firn profile (s).
!   grain_radius: Grain radius along the firn profile (m).
!
! Result:
!   d_density: Density change along the firn profile (kg m**-3).
! ----------------------------------------------------------------------

    ! doing nothing
    call tfm_essentials_do_nothing(nz, age)
    call tfm_essentials_do_nothing(nz, grain_radius)

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
! ----------------------------------------------------------------------
! Function: tfm_density_greve2009
!
! Greve, R. and Blatter, H. Dynamics of Ice Sheets and Glaciers.
! Springer, Berlin, (2009).
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of variables "depth", "density", "temperature", "age",
!       "grain_radius", "d_density".
!   dt: Time step (s).
!   depth: Depth along the firn profile (m).
!   density: Density along the firn profile (kg m**-3).
!   temperature: Temeperature along the firn profile (K).
!   age: Age along the firn profile (s).
!   grain_radius: Grain radius along the firn profile (m).
!
! Result:
!   d_density: Density change along the firn profile (kg m**-3).
!-----------------------------------------------------------------------

    ! doing nothing
    call tfm_essentials_do_nothing(nz, age)
    call tfm_essentials_do_nothing(nz, grain_radius)

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
! ----------------------------------------------------------------------
! Function: tfm_density_zwinger2007
!
! Zwinger, T., Greve, R., Gagliardini, O., Shiraiwa, T., and Lyly, M. A
! full Stokes-flow thermo-mechanical model for firn and ice applied to
! the Gorshkov crater glacier, Kamchatka. Annals of Glaciology, 45,
! 29-37, (2007). https://doi.org/10.3189/172756407782282543
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of variables "depth", "density", "temperature", "age",
!       "grain_radius", "d_density".
!   dt: Time step (s).
!   depth: Depth along the firn profile (m).
!   density: Density along the firn profile (kg m**-3).
!   temperature: Temeperature along the firn profile (K).
!   age: Age along the firn profile (s).
!   grain_radius: Grain radius along the firn profile (m).
!
! Result:
!   d_density: Density change along the firn profile (kg m**-3).
!-----------------------------------------------------------------------

    ! doing nothing
    call tfm_essentials_do_nothing(nz, age)
    call tfm_essentials_do_nothing(nz, grain_radius)

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
! ----------------------------------------------------------------------
! Function: tfm_density_breant2017
!
! Breant, C., Martinerie, P., Orsi, A., Arnaud, L., and Landais, A.
! Modelling firn thickness evolution during the last deglaciation:
! constraints on sensitivity to temperature and impurities. Clim. Past,
! 13, 833-853, (2017). https://doi.org/10.5194/cp-13-833-2017
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of variables "depth", "density", "temperature", "age",
!       "grain_radius", "d_density".
!   dt: Time step (s).
!   depth: Depth along the firn profile (m).
!   density: Density along the firn profile (kg m**-3).
!   temperature: Temeperature along the firn profile (K).
!   age: Age along the firn profile (s).
!   grain_radius: Grain radius along the firn profile (m).
!
! Result:
!   d_density: Density change along the firn profile (kg m**-3).
!-----------------------------------------------------------------------

    call tfm_essentials_do_nothing(nz, age)
    call tfm_essentials_do_nothing(nz, grain_radius)

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
! ----------------------------------------------------------------------
! Function tfm_density_medley2020
!
! Medley, B., Neumann, T. A., Zwally, H. J., and Smith, B. E. Forty-year
! Simulations of Firn Processes over the Greenland and Antarctic Ice
! Sheets. The Cryosphere Discuss. [preprint], in review, (2020).
! https://doi.org/10.5194/tc-2020-266
!
! NOTE! This functions is based on the preprint version of the paper.
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of variables "depth", "density", "temperature", "age",
!       "grain_radius", "d_density".
!   dt: Time step (s).
!   depth: Depth along the firn profile (m).
!   density: Density along the firn profile (kg m**-3).
!   temperature: Temeperature along the firn profile (K).
!   age: Age along the firn profile (s).
!   grain_radius: Grain radius along the firn profile (m).
!
! Result:
!   d_density: Density change along the firn profile (kg m**-3).
! ----------------------------------------------------------------------

    call tfm_essentials_do_nothing(nz, grain_radius)

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
! ----------------------------------------------------------------------
! Function: tfm_density_herron1980
!
! Herron, M. M. and Langway, C. C. Firn Densification: An Empirical
! Model. Journal of Glaciology, 25 (93), 373-385, (1980).
! https://doi.org/10.3189/S0022143000015239
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of variables "depth", "density", "temperature", "age",
!       "grain_radius", "d_density".
!   dt: Time step (s).
!   depth: Depth along the firn profile (m).
!   density: Density along the firn profile (kg m**-3).
!   temperature: Temeperature along the firn profile (K).
!   age: Age along the firn profile (s).
!   grain_radius: Grain radius along the firn profile (m).
!
! Result:
!   d_density: Density change along the firn profile (kg m**-3).
! ----------------------------------------------------------------------

    call tfm_essentials_do_nothing(nz, grain_radius)

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


  function tfm_density_sintering(nz, dt, depth, density, temperature, &
    & age, grain_radius) result(d_density)
    implicit none

    integer, intent(in)                   :: nz
    real(prec), intent(in)                :: dt
    real(prec), dimension(nz), intent(in) :: depth
    real(prec), dimension(nz), intent(in) :: density
    real(prec), dimension(nz), intent(in) :: temperature
    real(prec), dimension(nz), intent(in) :: age
    real(prec), dimension(nz), intent(in) :: grain_radius
    real(prec), dimension(nz)             :: d_density

    integer                   :: n
    real(prec)                :: neck_radius
    real(prec)                :: pore_radius
    real(prec), dimension(nz) :: stress
    real(prec), dimension(nz) :: strain_rate

    real(prec), parameter :: AMPLITUDE_GBO = 4.25E-6_prec ! m
    real(prec), parameter :: MU = 0.7_prec


    call tfm_essentials_do_nothing(nz, age)

    ! computation of the current stress state
    call tfm_density_computeStress(nz, depth, density, stress)

    strain_rate(:) = 0.0_prec
    do n = 1, nz, 1

      strain_rate(n) = strain_rate(n) + tfm_density_CobleCreep( &
      &  temperature(n),                                        &
      &  grain_radius(n),                                       &
      &  stress(n)                                              &
      )
     
      if ( density(n) <= CRITICAL_DENSITY ) then

        neck_radius = (MU * grain_radius(n))

        strain_rate(n) = strain_rate(n) + tfm_density_GrainBoundarySliding( &
        &  temperature(n),                                                  &
        &  density(n),                                                      &
        &  grain_radius(n),                                                 &
        &  neck_radius,                                                     &
        &  AMPLITUDE_GBO,                                                   &
        &  stress(n)                                                        &
        )
      end if


      if ( density(n) <= CLOSEOFF_DENSITY ) then

        strain_rate(n) = strain_rate(n) + tfm_density_NabarroHerringCreep( &
        &  temperature(n),                                                 &
        &  density(n),                                                     &
        &  grain_radius(n),                                                &
        &  stress(n)                                                       &
        )

        strain_rate(n) = strain_rate(n) + tfm_density_DislocationCreep( &
        &  temperature(n),                                              &
        &  density(n),                                                  &
        &  stress(n)                                                    &
        )

      else if ( (density(n) > CLOSEOFF_DENSITY) .and. (density(n) < ICE_DENSITY) ) then
        
        pore_radius = (                                                     &
        &  ((1.0_prec - (density(n) / ICE_DENSITY))**(1.0_prec / 3.0_prec)) &
        &  * grain_radius(n)                                                &
        )
 
        strain_rate(n) = strain_rate(n) + tfm_density_NabarroHerringCreepMod( &
        &  temperature(n),                                                    &
        &  density(n),                                                        &
        &  grain_radius(n),                                                   &
        &  pore_radius,                                                       &
        &  stress(n)                                                          &
        )

        strain_rate(n) = strain_rate(n) + tfm_density_DislocationCreepMod( &
        &  temperature(n),                                                 &
        &  density(n),                                                     &
        &  stress(n)                                                       &
        )

      else if ( density(n) >= ICE_DENSITY ) then
        
        strain_rate(n) = strain_rate(n) + 0.0_prec

      else
        print *, 'Module: tfm_density'
        print *, 'Function: tfm_density_sintering'
        print *, ''
        print *, 'The density seems to show and irregular value.'
        print *, 'Stopping right here!'
        STOP
      end if
    end do

    ! density change
    d_density = -(dt * strain_rate * density)
  end function tfm_density_sintering


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
    real(prec)                :: driving_force
    real(prec)                :: neck_radius
    real(prec)                :: pore_radius

    real(prec), parameter :: AMPLITUDE_GBO = 4.0E-6_prec ! m
    real(prec), parameter :: MU = 0.7_prec
    real(prec), parameter :: AIR_PRESSURE = 101325.0 ! (Pa)
! ----------------------------------------------------------------------
! Function: tfm_density_arthern1998
!
! Arthern, R. J. and Wingham, D. J. The Natural Fluctuations of Firn
! Densification and Their Effect on the Geodetic Determination of Ice
! Sheet Mass Balance. Climatic Change, 40, 605-624, (1998).
! https://doi.org/10.1023/A:1005320713306
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of variables "depth", "density", "temperature", "age",
!       "grain_radius", "d_density".
!   dt: Time step (s).
!   depth: Depth along the firn profile (m).
!   density: Density along the firn profile (kg m**-3).
!   temperature: Temeperature along the firn profile (K).
!   age: Age along the firn profile (s).
!   grain_radius: Grain radius along the firn profile (m).
!
! Parameters:
!   AMPLITUDE_GBO: Amplitude of grain boundary obstructions (m), see
!     Arthern & Wingham (1998).
!   MU: Ratio of neck radius to grain radius.
!   AIR_PRESSURE: Normal air pressure for the calculation of the
!     Boyle-Mariotte law in the third stage of densification (Pa).
!
! Result:
!   d_density: Density change along the firn profile (kg m**-3).
! ----------------------------------------------------------------------

    call tfm_essentials_do_nothing(nz, age)

    ! computation of the current stress state
    call tfm_density_computeStress(nz, depth, density, stress)

    ! computation of the strain rate
    strain_rate(:) = 0.0_prec

    do n = 1, nz, 1
      
      if ( (density(n) > 0.0_prec) .and. (density(n) <= CRITICAL_DENSITY) ) then
        
        neck_radius = (MU * grain_radius(n))

        strain_rate(n) = strain_rate(n) + tfm_density_GrainBoundarySliding( &
        &  temperature(n),                                                  &
        &  density(n),                                                      &
        &  grain_radius(n),                                                 &
        &  neck_radius,                                                     &
        &  AMPLITUDE_GBO,                                                   &
        &  stress(n)                                                        &
        )

        print *, n, 'I', density(n), strain_rate(n), -dt * strain_rate(n) * density(n)

      else if ( (density(n) > CRITICAL_DENSITY) .and. (density(n) <= CLOSEOFF_DENSITY) ) then
        
        strain_rate(n) = strain_rate(n) + tfm_density_NabarroHerringCreep( &
        &  temperature(n),                                                 &
        &  density(n),                                                     &
        &  grain_radius(n),                                                &
        &  stress(n)                                                       &
        )

        strain_rate(n) = strain_rate(n) + tfm_density_CobleCreep( &
        &  temperature(n),                                        &
        &  grain_radius(n),                                       &
        &  stress(n)                                              &
        )

        strain_rate(n) = strain_rate(n) + tfm_density_DislocationCreep( &
        &  temperature(n),                                              &
        &  density(n),                                                  &
        &  stress(n)                                                    &
        )

        print *, n, 'II', density(n), strain_rate(n), -dt * strain_rate(n) * density(n)

      else if ( (density(n) > CLOSEOFF_DENSITY) .and. (density(n) < ICE_DENSITY)) then
        
        pore_radius = (                                                     &
        &  ((1.0_prec - (density(n) / ICE_DENSITY))**(1.0_prec / 3.0_prec)) &
        &  * grain_radius(n)                                                &
        )
        
        !driving_force = (                                      &
        !&  stress(n)                                           &
        !&  - (AIR_PRESSURE * (                                 &
        !&    (density(n) * (ICE_DENSITY - CLOSEOFF_DENSITY))   &
        !&    / (CLOSEOFF_DENSITY * (ICE_DENSITY - density(n))) &
        !&  ))                                                  &
        !)
        driving_force = stress(n)
         
        strain_rate(n) = strain_rate(n) + tfm_density_NabarroHerringCreepMod( &
        &  temperature(n),                                                    &
        &  density(n),                                                        &
        &  grain_radius(n),                                                   &
        &  pore_radius,                                                       &
        &  driving_force                                                      &
        )

        strain_rate(n) = strain_rate(n) + tfm_density_CobleCreep( &
        &  temperature(n),                                        &
        &  grain_radius(n),                                       &
        &  driving_force                                          &
        )

        strain_rate(n) = strain_rate(n) + tfm_density_DislocationCreepMod( &
        &  temperature(n),                                                 &
        &  density(n),                                                     &
        &  driving_force                                                   &
        )

        print *, n, 'III', density(n), strain_rate(n), -dt * strain_rate(n) * density(n)

      else if ( density(n) >= ICE_DENSITY ) then

        strain_rate(n) = strain_rate(n) + 0.0_prec

      else

        print *, 'Module: tfm_density'
        print *, 'Function: tfm_density_arthern1998'
        print *, ''
        print *, 'The density seems to show an irregular value.'
        print *, 'Stopping right here!'
        print *, density
        STOP

      end if
    end do

    STOP

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
! ----------------------------------------------------------------------
! Function: tfm_density_li2003
!
! Li, J., Zwally, H. J., Corneja, H., and Yi, D. Seasonal variation of
! snow-surface elevation in North Greenland as modeled and detected by
! satellite radar altimetry. Annals of Glaciology, 37, 223-238, (2003).
! https://doi.org/10.3189/172756403781815889
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of variables "depth", "density", "temperature", "age",
!       "grain_radius", "d_density".
!   dt: Time step (s).
!   depth: Depth along the firn profile (m).
!   density: Density along the firn profile (kg m**-3).
!   temperature: Temeperature along the firn profile (K).
!   age: Age along the firn profile (s).
!   grain_radius: Grain radius along the firn profile (m).
!
! Result:
!   d_density: Density change along the firn profile (kg m**-3).
! ----------------------------------------------------------------------

    call tfm_essentials_do_nothing(nz, grain_radius)

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
! ----------------------------------------------------------------------
! Function: tfm_density_helsen2008
!
! Helsen, M. M., van den Broeke, M. R., van den Wal, R. S. W.,
! van de Berg, W. J., van Meijgaard, E., Davis, C. H., Li, Y., and
! Goodwin, I. Elevation Changes in Antarctica Mainly Determined by
! Accumulation Variability. Science, 320 (5883), 1626-1629, (2008).
! https://doi.org/10.1126/science.1153894
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of variables "depth", "density", "temperature", "age",
!       "grain_radius", "d_density".
!   dt: Time step (s).
!   depth: Depth along the firn profile (m).
!   density: Density along the firn profile (kg m**-3).
!   temperature: Temeperature along the firn profile (K).
!   age: Age along the firn profile (s).
!   grain_radius: Grain radius along the firn profile (m).
!
! Result:
!   d_density: Density change along the firn profile (kg m**-3).
! ----------------------------------------------------------------------

    call tfm_essentials_do_nothing(nz, grain_radius)

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
! ----------------------------------------------------------------------
! Function: tfm_density_arthern2010
!
! Arthern, R. J., Vaughan, D. G., Rankin, A. M., Mulvaney, R., and
! Thomas, E. R. In situ measurements of Antarctic snow compaction
! compared with predictions of models. Journal of Geophysical Research:
! Earth Surface, 115 (F3), (2010). https://doi.org/10.1029/2009JF001306
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of variables "depth", "density", "temperature", "age",
!       "grain_radius", "d_density".
!   dt: Time step (s).
!   depth: Depth along the firn profile (m).
!   density: Density along the firn profile (kg m**-3).
!   temperature: Temeperature along the firn profile (K).
!   age: Age along the firn profile (s).
!   grain_radius: Grain radius along the firn profile (m).
!
! Result:
!   d_density: Density change along the firn profile (kg m**-3).
! ----------------------------------------------------------------------

    call tfm_essentials_do_nothing(nz, grain_radius)

    ! 10 m temperature
    do n = nz, 1, -1
      if ( (depth(n) - depth(nz)) <= -10.0_prec ) EXIT
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
! ----------------------------------------------------------------------
! Function: tfm_density_ligtenberg2011
!
! Ligtenberg, S. R. M., Helsen, M. M. and van den Broeke, M. R. An
! improved semi-empirical model for the densification of Antarctic firn.
! The Cryosphere, 5, 809-819, (2011).
! https://doi.org/10.5194/tc-5-809-2011
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of variables "depth", "density", "temperature", "age",
!       "grain_radius", "d_density".
!   dt: Time step (s).
!   depth: Depth along the firn profile (m).
!   density: Density along the firn profile (kg m**-3).
!   temperature: Temeperature along the firn profile (K).
!   age: Age along the firn profile (s).
!   grain_radius: Grain radius along the firn profile (m).
!
! Result:
!   d_density: Density change along the firn profile (kg m**-3).
! ----------------------------------------------------------------------

    d_density = tfm_density_arthern2010(nz, dt, depth, density, &
      & temperature, age, grain_radius)

    ! boundary between first and seconds stage
    do n = nz, 1, -1
      if ( density(n) >= 550.0 ) EXIT
    end do

    ! change accumulation to (kg a-1 m-2)
    call tfm_essentials_mean_acc(nz, depth, density, age, mean_acc)
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
! ----------------------------------------------------------------------
! Function: tfm_density_simonsen2013
!
! Simonsen, S. B., Stenseng, L., Adalgeisdottir, G., Fausto, R. S.,
! Hvidberg, C. S., and Lucas-Picher, P. Assessing a multilayered dynamic
! firn-compaction model for Greenland with ASIRAS radar measurements.
! Journal of Glaciology, 59 (215), 545-558, (2013).
! https://doi.org/10.3189/2013JoG12J158
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of variables "depth", "density", "temperature", "age",
!       "grain_radius", "d_density".
!   dt: Time step (s).
!   depth: Depth along the firn profile (m).
!   density: Density along the firn profile (kg m**-3).
!   temperature: Temeperature along the firn profile (K).
!   age: Age along the firn profile (s).
!   grain_radius: Grain radius along the firn profile (m).
!
! Result:
!   d_density: Density change along the firn profile (kg m**-3).
! ----------------------------------------------------------------------

    d_density = tfm_density_arthern2010(nz, dt, depth, density, &
      & temperature, age, grain_radius)

    ! boundary between first and seconds stage
    do n = nz, 1, -1
      if ( density(n) >= 550.0 ) EXIT
    end do

    ! change accumulation to (kg a-1 m-2)
    call tfm_essentials_mean_acc(nz, depth, density, age, mean_acc)
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
! ---------------------------------------------------------------------
! Subroutine: tfm_density_herron1980_analytical
!
! Herron, M. M. and Langway, C. C. Firn Densification: An Empirical
! Model. Journal of Glaciology, 25 (93), 373-385, (1980).
! https://doi.org/10.3189/S0022143000015239
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of variables "depth" and "density".
!   temperature: Temperature (K).
!   accumulation: Accumulation rate (m weq. s**-1).
!   surface_density: Surface density (kg m**-3).
!   depth: Depth along the firn profile (m).
!   density - on input: Variable to store the result.
!
! Result:
!   density - on output: Density (kg m**-3).
! ---------------------------------------------------------------------

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
