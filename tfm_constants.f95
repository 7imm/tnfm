module tfm_constants
  use settings
  implicit none

  ! acceleration due to gravity (kg / m s**2)
  real(prec), parameter :: ACC_GRAVITY = 9.81

  ! gas constant (J / K mol)
  real(prec), parameter :: GAS_CONST = 8.31446261815324

  ! ice density (kg / m**3)
  real(prec), parameter :: ICE_DENSITY = 917.0

  ! water density (kg / m**3)
  real(prec), parameter :: WATER_DENSITY = 1000.0

  ! seconds per year (s / yr)
  real(prec), parameter :: SECONDS_YEAR = (3600.0 * 24.0 * 365.0)

  ! seconds per day (s / d)
  real(prec), parameter :: SECONDS_DAY = (3600.0 * 24.0)

  ! specific heat capacity of ice (J / kg K) Reijmer et al. 2012
  real(prec), parameter :: SPECIFIC_HEAT_ICE = 2050.0

  ! latent heat of ice (J / kg) Reijmer et al. 2012
  real(prec), parameter :: LATENT_HEAT = 334000.0
end module tfm_constants
