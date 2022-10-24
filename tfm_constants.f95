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
  real(prec), parameter :: BOLTZMANN = 1.380649e-23_prec

  ! ice flow exponent
  real(prec), parameter :: ICE_N = 3.0_prec

  ! bubble close off density
  real(prec), parameter :: CLOSEOFF_DENSITY = 834.0_prec

  ! tempreature offset to convert K and C
  real(prec), parameter :: TEMP_OFFSET = 273.15_prec

  ! pi
  real(prec), parameter :: PI = 4.0d0 * datan(1.d0)
end module tfm_constants
