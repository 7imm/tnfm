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
!  tfm_grain_zwally2002: Zwally & Li (2002)
!  tfm_grain_brun1989: Brun (1989)
!  tfm_grain_tusima1978: Tusima (1978)
!  tfm_grain_katsushima2009: Katsushima et al. (2009)
! ----------------------------------------------------------------------

  contains


  subroutine tfm_grain_Area2RadiusChangeExplicit(nz, dt, grain_radius, &
    & area_growth_rate, d_grain_radius)
    implicit none

    integer, intent(in)                      :: nz
    real(prec), intent(in)                   :: dt
    real(prec), dimension(nz), intent(in)    :: grain_radius
    real(prec), dimension(nz), intent(in)    :: area_growth_rate
    real(prec), dimension(nz), intent(inout) :: d_grain_radius

    d_grain_radius = (                &
    &  ((2.0_prec * PI)**(-1.0_prec)) &
    &  * (grain_radius**(-1.0_prec))  &
    &  * area_growth_rate             &
    &  * dt                           &
    )
  end subroutine tfm_grain_Area2RadiusChangeExplicit


  subroutine tfm_grain_Area2RadiusChangeImplicit(nz, dt, grain_radius, &
    & area_growth_rate, d_grain_radius)
    implicit none

    integer, intent(in)                      :: nz
    real(prec), intent(in)                   :: dt
    real(prec), dimension(nz), intent(in)    :: grain_radius
    real(prec), dimension(nz), intent(in)    :: area_growth_rate
    real(prec), dimension(nz), intent(inout) ::d_grain_radius

    d_grain_radius = ((                   &
    &  (grain_radius / 2.0_prec)          &
    &  + ((                               &
    &    ((grain_radius**2.0_prec) / 4.0) &
    &    + (                              &
    &      ((2.0_prec * PI)**(-1.0_prec)) &
    &      * area_growth_rate             &
    &      * dt                           &
    &    )                                &
    &  )**0.5_prec)                       &
    ) - grain_radius)
  end subroutine tfm_grain_Area2RadiusChangeImplicit


  subroutine tfm_grain_Volume2RadiusChangeExplicit(nz, dt, grain_radius, &
    & volume_growth_rate, d_grain_radius)
    implicit none

    integer, intent(in)                      :: nz
    real(prec), intent(in)                   :: dt
    real(prec), dimension(nz), intent(in)    :: grain_radius
    real(prec), dimension(nz), intent(in)    :: volume_growth_rate
    real(prec), dimension(nz), intent(inout) :: d_grain_radius

    d_grain_radius = (                &
    &  ((4.0_prec * PI)**(-1.0_prec)) &
    &  * (grain_radius**(-2.0_prec))  &
    &  * volume_growth_rate           &
    &  * dt                           &
    )
  end subroutine tfm_grain_Volume2RadiusChangeExplicit


  subroutine tfm_grain_Volume2RadiusChangeImplicit(nz, dt, grain_radius, &
    & volume_growth_rate, d_grain_radius)
    implicit none

    integer, intent(in)                      :: nz
    real(prec), intent(in)                   :: dt
    real(prec), dimension(nz), intent(in)    :: grain_radius
    real(prec), dimension(nz), intent(in)    :: volume_growth_rate
    real(prec), dimension(nz), intent(inout) :: d_grain_radius

    real(prec), dimension(nz) :: p
    real(prec), dimension(nz) :: q
    real(prec), dimension(nz) :: delta
    real(prec), dimension(nz) :: u
    real(prec), dimension(nz) :: v

    p = ((-1.0_prec / 3.0_prec) * (-grain_radius**2.0_prec))
    q = (                                                          &
    &  ((2.0_prec / 27.0_prec) * (-grain_radius**3.0))             &
    &  + ((-1.0_prec / (4.0_prec * PI)) * volume_growth_rate * dt) &
    )

    delta = (((q**2.0_prec) / 4.0_prec) + ((p**3.0_prec) / 27.0_prec))

    u = (((-q / 2.0_prec) + (delta**0.5_prec))**(1.0_prec / 3.0_prec))
    v = (((-q / 2.0_prec) - (delta**0.5_prec))**(1.0_prec / 3.0_prec))

    d_grain_radius = ((grain_radius / 3.0_prec) + u + v - grain_radius)
  end subroutine tfm_grain_Volume2RadiusChangeImplicit


  subroutine tfm_grain_Area2RadiusChange(nz, dt, grain_radius, &
    & area_growth_rate, d_grain_radius, solving_method)
    implicit none

    integer, intent(in)                      :: nz
    real(prec), intent(in)                   :: dt
    real(prec), dimension(nz), intent(in)    :: grain_radius
    real(prec), dimension(nz), intent(in)    :: area_growth_rate
    real(prec), dimension(nz), intent(inout) :: d_grain_radius
    character(len=*), intent(in)             :: solving_method

    if ( solving_method == 'explicit' ) then
      call tfm_grain_Area2RadiusChangeExplicit(nz, dt, grain_radius, &
        & area_growth_rate, d_grain_radius)

    else if ( solving_method == 'implicit' ) then
      call tfm_grain_Area2RadiusChangeImplicit(nz, dt, grain_radius, &
        & area_growth_rate, d_grain_radius)

    else
      print *, '*******************************************************'
      print *, '* Module: tfm_grain                                   *'
      print *, '* Subroutine: tfm_grain_Area2RadiusChange             *'
      print *, '*                                                     *'
      print *, '* Argument "solving_method" has eiter to be           *'
      print *, '* "explicit" or "implicit". This seems not to be the  *'
      print *, '* case.                                               *'
      print *, '* Stopping right here!                                *'
      print *, '*******************************************************'
      STOP
    end if
  end subroutine tfm_grain_Area2RadiusChange


  subroutine tfm_grain_Volume2RadiusChange(nz, dt, grain_radius, &
    & volume_growth_rate, d_grain_radius, solving_method)
    implicit none
    
    integer, intent(in)                      :: nz
    real(prec), intent(in)                   :: dt
    real(prec), dimension(nz), intent(in)    :: grain_radius
    real(prec), dimension(nz), intent(in)    :: volume_growth_rate
    real(prec), dimension(nz), intent(inout) :: d_grain_radius
    character(len=*), intent(in)             :: solving_method

    if ( solving_method == 'explicit' ) then
      call tfm_grain_Volume2RadiusChangeExplicit(nz, dt, grain_radius, &
        & volume_growth_rate, d_grain_radius)

    else if ( solving_method == 'implicit' ) then
      call tfm_grain_Volume2RadiusChangeImplicit(nz, dt, grain_radius, &
        & volume_growth_rate, d_grain_radius)
      
    else
      print *, '*******************************************************'
      print *, '* Module: tfm_grain                                   *'
      print *, '* Subroutine: tfm_grain_Volume2RadiusChange           *'
      print *, '*                                                     *'
      print *, '* Argument "solving_method" has eiter to be           *'
      print *, '* "explicit" or "implicit". This seems not to be the  *'
      print *, '* case.                                               *'
      print *, '* Stopping right here!                                *'
      print *, '*******************************************************'
      STOP
    end if
  end subroutine tfm_grain_Volume2RadiusChange


  function tfm_grain_arthern2010(nz, dt, temperature, density, &
    & liquid_water, grain_radius) result(d_grain_radius)
    implicit none

    ! Parameters from Arthern & Wingham (2010) following Paterson (1994)
    real(prec), parameter :: PRE_FACTOR = 1.3e-7_prec ! m**2 s**-1
    real(prec), parameter :: ACTIVATION_ENERGY = 42400.0_prec !J mol**-1

    integer, intent(in)                   :: nz
    real(prec), intent(in)                :: dt
    real(prec), dimension(nz), intent(in) :: temperature
    real(prec), dimension(nz), intent(in) :: density
    real(prec), dimension(nz), intent(in) :: liquid_water
    real(prec), dimension(nz), intent(in) :: grain_radius

    real(prec), dimension(nz) :: d_grain_radius
    real(prec), dimension(nz) :: area_growth_rate
! ----------------------------------------------------------------------
! Function: tfm_grain_arthern2010
!
! Arthern, R. J., Vaughan, D. G., Rankin, A. M., Mulvaney, R., and
! Thomas, E. R. In situ measurements of Antarctic snow compaction
! compared with predicitions of models. Journal of Geophysical Research:
! Earth Surface, 115 (F3), (2010). https://doi.org/10.1029/2009JF001306
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of the variables "temperature", "density", and
!       "d_grain_radius".
!   dt: Time step (s).
!   temperature: Temperature along the firn profile (K).
!   density: Density along the firn profile (kg / m**3).
!   liquid_water: Volumetric liquid water content along the profile (1).
!   grain_radius: Grain radius alon ght eprofile (m).
!
! Result:
!   d_grain_radius: Grain radius change along the firn profile (m).
! ----------------------------------------------------------------------

    call tfm_essentials_do_nothing(nz, density)
    call tfm_essentials_do_nothing(nz, liquid_water)
    
    ! growth rate (m**2 s**-1)
    area_growth_rate = (                                     &
    &  PRE_FACTOR                                            &
    &  * exp(-ACTIVATION_ENERGY / (GAS_CONST * temperature)) &
    )

    call tfm_grain_Area2RadiusChange(nz, dt, grain_radius, &
      & area_growth_rate, d_grain_radius, 'implicit')
  end function tfm_grain_arthern2010


  function tfm_grain_zwally2002(nz, dt, temperature, density, &
    & liquid_water, grain_radius) result(d_grain_radius)
    implicit none

    integer, intent(in)                   :: nz
    real(prec), intent(in)                :: dt
    real(prec), dimension(nz), intent(in) :: temperature
    real(prec), dimension(nz), intent(in) :: density
    real(prec), dimension(nz), intent(in) :: liquid_water
    real(prec), dimension(nz), intent(in) :: grain_radius

    real(prec), dimension(nz) :: area_growth_rate
    real(prec), dimension(nz) :: d_grain_radius
! ----------------------------------------------------------------------
! Function: tfm_grain_li2002
!
! Grain growth following Zwally & Li (2002) as described in
! Li & Zwally (2011).
!
! Zwally, H. J. and Li, J. Seasonal and interannual variations of firn
! densification and ice-sheet elevation at the Greenland summit.
! Journal of Glaciology, 48 (171), (2002).
! https://doi.org/10.3189/172756502781831403
!
! Li, J. and Zwally, H. J. Modeling of firn compaction for estimating
! ice-sheet mass change from observed ice-sheet elevation change. Annals
! of Glaciology, 52 (59), (2011).
! https://doi.org/10.3189/172756411799096321
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of the variables "temperature", "density", and
!       "d_grain_radius".
!   dt: Time step (s).
!   temperature: Temperature along the firn profile (K).
!   density: Density along the firn profile (kg / m**3).
!   liquid_water: Volumetric liquid water content along the profile (1).
!   grain_radius: Grain radius alon ght eprofile (m).
!
! Result:
!   d_grain_radius: Grain radius change along the firn profile (m).
! ----------------------------------------------------------------------
    
    call tfm_essentials_do_nothing(nz, density)
    call tfm_essentials_do_nothing(nz, liquid_water)

    ! growth rate (mm**2 a**-1)
    area_growth_rate = 8.36_prec * ((273.2_prec - temperature)**(-2.061_prec))
    area_growth_rate = area_growth_rate * 1.0E-6_prec / SECONDS_YEAR ! (m**2 s**-1)

    call tfm_grain_Area2RadiusChange(nz, dt, grain_radius, &
      & area_growth_rate, d_grain_radius, 'implicit')
  end function tfm_grain_zwally2002
  

  function tfm_grain_brun1989(nz, dt, temperature, density, &
    & liquid_water, grain_radius) result(d_grain_radius)
    implicit none

    integer, intent(in)                   :: nz
    real(prec), intent(in)                :: dt
    real(prec), dimension(nz), intent(in) :: temperature
    real(prec), dimension(nz), intent(in) :: density
    real(prec), dimension(nz), intent(in) :: liquid_water
    real(prec), dimension(nz), intent(in) :: grain_radius

    real(prec), dimension(nz) :: d_grain_radius
    real(prec), dimension(nz) :: gravimetric_wc
    real(prec), dimension(nz) :: volume_growth_rate
! ----------------------------------------------------------------------
! Fucntion: tfm_grain_brun1989
!
! Grain growth following Brun 1989.
!
! Brun, E. (1989). Investigation On Wet-Snow Metamorphism in Respect of
! Liquid-Water Content. Annals of Glaciology, 13, pp. 22-26,
! https://doi.org/10.3189/S0260305500007576
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of the variables "temperature", "density", and
!       "d_grain_radius".
!   dt: Time step (s).
!   temperature: Temperature along the firn profile (K).
!   density: Density along the firn profile (kg / m**3).
!   liquid_water: Volumetric liquid water content along the profile (1).
!   grain_radius: Grain radius alon ght eprofile (m).
!
! Result:
!   d_grain_radius: Grain radius change along the firn profile (m).
! ----------------------------------------------------------------------
    
    call tfm_essentials_do_nothing(nz, temperature)

    gravimetric_wc = (100.0_prec * (                &
    &  (liquid_water * WATER_DENSITY)               &
    &  / (density + (liquid_water * WATER_DENSITY)) &
    ))
    
    ! growth rate (mm**3 s**-1)
    volume_growth_rate = (1.28E-8_prec + (4.22E-10_prec * (gravimetric_wc**3.0_prec)))
    volume_growth_rate = (volume_growth_rate * 1.0E-9_prec) ! (m**3 s**-1)

    call tfm_grain_Volume2RadiusChange(nz, dt, grain_radius, &
      & volume_growth_rate, d_grain_radius, 'explicit')
  end function tfm_grain_brun1989


  function tfm_grain_tusima1978(nz, dt, temperature, density, &
    & liquid_water, grain_radius) result(d_grain_radius)
    implicit none

    integer, intent(in)    :: nz
    real(prec), intent(in) :: dt
    real(prec), dimension(nz), intent(in) :: temperature
    real(prec), dimension(nz), intent(in) :: density
    real(prec), dimension(nz), intent(in) :: liquid_water
    real(prec), dimension(nz), intent(in) :: grain_radius

    real(prec), dimension(nz) :: d_grain_radius
    real(prec), dimension(nz) :: radius_growth_rate 
! ----------------------------------------------------------------------
! Function tfm_grain_tusima1978
!
! Grain growth following Tusima (1978) as described by Katsushima et al. 
! (2009). The original paper is written in Japanese.
!
! Tusima, K. (1978). Grain coarsening of ice particles immersed in pure
! water. Seppyo, 40 (4), pp. 155-165.
!
! Katsushima, T., Kamakura, T., and Takeuchi, Y. (2009). A multiple
! layer model including a parametrization of vertical water channel
! process in snowpack. Cold Regions Science and Technology, 59, pp.
! 143-151, https://doi.org/10.1016/j.coldregions.2009.09.002
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of the variables "temperature", "density", and
!       "d_grain_radius".
!   dt: Time step (s).
!   temperature: Temperature along the firn profile (K).
!   density: Density along the firn profile (kg / m**3).
!   liquid_water: Volumetric liquid water content along the profile (1).
!   grain_radius: Grain radius alon ght eprofile (m).
!
! Result:
!   d_grain_radius: Grain radius change along the firn profile (m).
! ----------------------------------------------------------------------
    
    call tfm_essentials_do_nothing(nz, temperature)
    call tfm_essentials_do_nothing(nz, density)
    call tfm_essentials_do_nothing(nz, liquid_water)

    ! growth rate (mm s**-1)
    radius_growth_rate = (                                     &
    &  (2.5E-4_prec / ((grain_radius * 1.0E3_prec)**2.0_prec)) &
    &  / 3600.0_prec                                           &
    )
    radius_growth_rate = (radius_growth_rate * 1.0E-3_prec) ! (m s**-1)

    d_grain_radius = (dt * radius_growth_rate)
  end function tfm_grain_tusima1978


  function tfm_grain_katsushima2009(nz, dt, temperature, density, &
    & liquid_water, grain_radius) result(d_grain_radius)
    implicit none

    integer, intent(in)                   :: nz
    real(prec), intent(in)                :: dt
    real(prec), dimension(nz), intent(in) :: temperature
    real(prec), dimension(nz), intent(in) :: density
    real(prec), dimension(nz), intent(in) :: liquid_water
    real(prec), dimension(nz), intent(in) :: grain_radius

    real(prec), dimension(nz) :: gravimetric_wc
    real(prec), dimension(nz) :: d_grain_radius_brun
    real(prec), dimension(nz) :: d_grain_radius_tusima
    real(prec), dimension(nz) :: d_grain_radius
! ----------------------------------------------------------------------
! Function tfm_grain_Katushima2009
!
! Grain growth following Katsushima et al. (2009) combining the
! parametrizations of Brun (1989) and Tusima (1978) depending on the
! liquid water content of the firn.
!
! Katsushima, T., Kamakura, T., and Takeuchi, Y. (2009). A multiple
! layer model including a parametrization of vertical water channel
! process in snowpack. Cold Regions Science and Technology, 59, pp.
! 143-151, https://doi.org/10.1016/j.coldregions.2009.09.002
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of the variables "temperature", "density", and
!       "d_grain_radius".
!   dt: Time step (s).
!   temperature: Temperature along the firn profile (K).
!   density: Density along the firn profile (kg / m**3).
!   liquid_water: Volumetric liquid water content along the profile (1).
!   grain_radius: Grain radius alon ght eprofile (m).
!
! Result:
!   d_grain_radius: Grain radius change along the firn profile (m).
! ----------------------------------------------------------------------

    call tfm_essentials_do_nothing(nz, temperature)

    gravimetric_wc = (100.0_prec * (                &
    &  (liquid_water * WATER_DENSITY)               &
    &  / (density + (liquid_water * WATER_DENSITY)) &
    ))
 
    d_grain_radius_brun = tfm_grain_brun1989(nz, dt, temperature, density, &
      & liquid_water, grain_radius)

    d_grain_radius_tusima = tfm_grain_tusima1978(nz, dt, temperature, density, &
      & liquid_water, grain_radius)

    where ( gravimetric_wc <= 10.0_prec )
      d_grain_radius = d_grain_radius_brun
    else where
      d_grain_radius = min(d_grain_radius_brun, d_grain_radius_tusima)
    end where
  end function tfm_grain_katsushima2009
end module tfm_grain
