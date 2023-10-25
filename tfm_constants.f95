module tfm_constants
  implicit none

  ! precision
  integer, parameter :: prec = selected_real_kind(8)

  ! acceleration due to gravity (kg / m s**2)
  real(prec), parameter :: ACC_GRAVITY = 9.81_prec

  ! gas constant (J / K mol)
  real(prec), parameter :: GAS_CONST = 8.31446261815324_prec

  ! ice density (kg / m**3)
  real(prec), parameter :: ICE_DENSITY = 917.0_prec

  ! water density (kg / m**3)
  real(prec), parameter :: WATER_DENSITY = 1000.0_prec

  ! seconds per year (s / yr)
  real(prec), parameter :: SECONDS_YEAR = (3600.0_prec * 24.0_prec * 365.0_prec)

  ! seconds per day (s / d)
  real(prec), parameter :: SECONDS_DAY = (3600.0_prec * 24.0_prec)

  ! specific heat capacity of ice (J / kg K) Reijmer et al. 2012
  real(prec), parameter :: SPECIFIC_HEAT_ICE = 2050.0_prec

  ! latent heat of ice (J / kg) Reijmer et al. 2012
  real(prec), parameter :: LATENT_HEAT = 334000.0_prec

  ! melt tempreature (K)
  real(prec), parameter :: MELT_TEMP = 273.15_prec

  ! Boltzmann constant
  real(prec), parameter :: BOLTZMANN = 1.380649E-23_prec

  ! ice flow exponent
  real(prec), parameter :: ICE_N = 3.0_prec

  ! critical density
  real(prec), parameter :: CRITICAL_DENSITY = 550.0_prec

  ! bubble close off density
  real(prec), parameter :: CLOSEOFF_DENSITY = 834.0_prec

  ! tempreature offset to convert K and C
  real(prec), parameter :: TEMP_OFFSET = 273.15_prec

  ! pi
  real(prec), parameter :: PI = 4.0d0 * datan(1.d0)

  ! Volume of a H2O molecule (Maeno & Ebinuma, 1983)
  real(prec), parameter :: VOLUME_H2O = 3.27E-29_prec

  ! Burgers vector of ice (Maeno & Ebinuma, 1983)
  real(prec), parameter :: BURGERS_VECTOR = 0.45E-9_prec ! m

  ! pre factor and activation energy for lattice diffusion in ice (Maeno & Ebinnuma, 1983)
  real(prec), parameter :: PRE_FACTOR_LATTICE_DIFFUSION = 0.03_prec
  real(prec), parameter :: ACTIVATION_ENERGY_LATTICE_DIFFUSION = 66200.0_prec ! J mol**-1

  ! pre factor and activatio energy for boundary diffusion in ice (Maeno & Ebinuma, 1983)
  real(prec), parameter :: PRE_FACTOR_BOUNDARY_DIFFUSION = 0.03_prec
  real(prec), parameter :: ACTIVATION_ENERGY_BOUNDARY_DIFFUSION = 44100.0_prec ! J mol**-1

  ! pre factor and activation energy for dislocation creep in ice (Greve & Blatter, 2009)
  real(prec), parameter :: PRE_FACTOR_DISLOCATION_CREEP_LOW = 3.985E-13_prec
  real(prec), parameter :: ACTIVATION_ENERGY_DISLOCATION_CREEP_LOW = 60000.0_prec
  real(prec), parameter :: PRE_FACTOR_DISLOCATION_CREEP_HIGH = 1.916E3_prec
  real(prec), parameter :: ACTIVATION_ENERGY_DISLOCATION_CREEP_HIGH = 139000.0_prec
end module tfm_constants
