module tfm_grain
  use tfm_essentials 
  use tfm_constants
  implicit none
! ----------------------------------------------------------------------
! Module: tfm_grain
!
! Dependencies: tfm_essentials, tfm_constants
!
! Functions:
!  tfm_grain_arthern2010: Arthern et al. (2010)
!
! Subroutines:
!  tfm_grain_gowType: Gow type grain growth.
! ----------------------------------------------------------------------

  contains


  subroutine tfm_grain_gowType(nz, dt, pre_factor, activation_energy, &
    & temperature, d_grain_radius)
    implicit none

    integer, intent(in)                      :: nz
    real(prec), intent(in)                   :: dt
    real(prec), intent(in)                   :: pre_factor
    real(prec), intent(in)                   :: activation_energy
    real(prec), dimension(nz), intent(in)    :: temperature
    real(prec), dimension(nz), intent(inout) :: d_grain_radius
! ----------------------------------------------------------------------
! Subroutine: tfm_grain_gowType
! 
! Generic function following the temperature dependent grain growth
! approach using an Arrhenius equation. The grain growth equation of Gow
! is formulated in a differential form (see Schultz et al. 2022).
!
! Gow, A. J. On the Rates of Growth of Grains and Crystals in South
! Polar Firn. Journal of Glaciology, 8 (53), 241-252, (1969).
! https://doi.org/10.3189/S0022143000031233
!
! Schultz, T., Müller, R., Gross, D., and Humbert, A. On the
! contribution of grain boundary sliding type creep to firn
! densification - an assessment using an optimization approach. The
! Cryosphere, 16, 143-158, (2022).https://doi.org/10.5194/tc-16-143-2022
!
! Author: Timm Schultz
! 
! Arguments:
!   nz: Dimension of variables "temperature", "grain_radius"
!     and "d_grain_radius".
!   dt: Time step (s).
!   temperature: Temperature along the firn profile (K).
!   pre_factor: Pre factor of the Arrhenius equation (m**2 s**-1).
!   activation_energy: Activation energy for the Arrhenius
!     equation (J mol**-1).
!   d_grain_radius (on input): Dummy variable for the change of the
!     grain radius.
!
! Result:
!  d_grain_radius (on output): Change of grain radius along the firn
!    profile (m).
! ----------------------------------------------------------------------
    
    d_grain_radius = (                                                   &
      & pre_factor * exp((-activation_energy)/(GAS_CONST * temperature)) &
    )
    d_grain_radius = dt * (d_grain_radius**0.5)
  end subroutine tfm_grain_gowType
  

  function tfm_grain_arthern2010(nz, dt, depth, density, temperature) &
    & result(d_grain_radius)
    implicit none

    ! Parameters from Arthern & Wingham (2010) following Paterson (1994)
    real(prec), parameter :: PRE_FACTOR = 1.3e-7_prec ! m**2 s**-1
    real(prec), parameter :: ACTIVATION_ENERGY = 42400.0_prec !J mol**-1

    integer, intent(in)                   :: nz
    real(prec), intent(in)                :: dt
    real(prec), dimension(nz), intent(in) :: depth
    real(prec), dimension(nz), intent(in) :: density
    real(prec), dimension(nz), intent(in) :: temperature

    real(prec), dimension(nz) :: d_grain_radius
! ----------------------------------------------------------------------
! Function: tfm_density_arthern2010
!
! Arthern, R. J., Vaughan, D. G., Rankin, A. M., Mulvaney, R., and
! Thomas, E. R. In situ measurements of Antarctic snow compaction
! compared with predicitions of models. Journal of Geophysical Research:
! Earth Surface, 115 (F3), (2010). https://doi.org/10.1029/2009JF001306
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of vairbales "temperature" and "grain radius".
!   dt: Time step (s).
!   temperature: Temperature along the firn profile (K).
!   grain_radius: Grain radius along the firn profile (m).
!
! Result:
!   d_grain_radius: Grain radius change along the firn profile (m).
! ----------------------------------------------------------------------
    
    call tfm_essentials_do_nothing(nz, depth)
    call tfm_essentials_do_nothing(nz, density)
    
    call tfm_grain_gowType( &
    &  nz,                  &
    &  dt,                  &
    &  PRE_FACtOR,          &
    &  ACTIVATION_ENERGY,   &
    &  temperature,         &
    &  d_grain_radius       &
    )
  end function tfm_grain_arthern2010


  function tfm_grain_li2002(nz, dt, depth, density, temperature) &
    & result(d_grain_radius)
    implicit none

    integer, intent(in)                   :: nz
    real(prec), intent(in)                :: dt
    real(prec), dimension(nz), intent(in) :: depth
    real(prec), dimension(nz), intent(in) :: density
    real(prec), dimension(nz), intent(in) :: temperature

    real(prec), dimension(nz) :: d_grain_radius

  end function tfm_grain_li2002
end module tfm_grain
