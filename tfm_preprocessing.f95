module tfm_preprocessing
  use tfm_essentials 
  use tfm_constants
  implicit none

  contains

  subroutine tfm_pre_pdd(nt, dt, temperature, parameters, melt)
    implicit none

    integer, intent(in)                   :: nt
    real(prec), intent(in)                :: dt
    real(prec), dimension(nt), intent(in) :: temperature
    character(len=*), intent(in)          :: parameters

    real(prec), dimension(nt), intent(inout) :: melt

    integer    :: n
    real(prec) :: temp0, ddf

    ! parameters
    if ( parameters == 'medley2020_greenland' ) then
      call tfm_pre_pdd_medley2020_greenland(ddf, temp0)

    else if ( parameters == 'ruckamp2019_greenland_snow' ) then
      call tfm_pre_pdd_ruckamp2019_greenland_snow(ddf, temp0)

    else if ( parameters == 'ruckamp2019_greenland_ice' ) then
      call tfm_pre_pdd_ruckamp2019_greenland_ice(ddf, temp0)

    else
      print *, 'no parameter set for PDD scheme given'
      STOP

    end if

    ! PDD scheme
    melt = 0.0
    do n = 1, nt, 1
      if ( temperature(n) > temp0 ) then
        melt(n) = (dt * ddf)
      end if
    end do
  end subroutine tfm_pre_pdd

  
  subroutine tfm_pre_pdd_medley2020_greenland(ddf, temp0)
    implicit none

    real(prec), intent(inout) :: ddf
    real(prec), intent(inout) :: temp0

    ddf = 54.0 / WATER_DENSITY / SECONDS_YEAR
    temp0 = 269.0
  end subroutine tfm_pre_pdd_medley2020_greenland

  
  subroutine tfm_pre_pdd_ruckamp2019_greenland_snow(ddf, temp0)
    implicit none

    real(prec), intent(inout) :: ddf
    real(prec), intent(inout) :: temp0

    ddf = 3.0e-3 / SECONDS_DAY
    temp0 = 273.16
  end subroutine tfm_pre_pdd_ruckamp2019_greenland_snow


  subroutine tfm_pre_pdd_ruckamp2019_greenland_ice(ddf, temp0)
    implicit none

    real(prec), intent(inout) :: ddf
    real(prec), intent(inout) :: temp0

    ddf = 8.0e-3 / SECONDS_DAY
    temp0 = 273.16
  end subroutine tfm_pre_pdd_ruckamp2019_greenland_ice


  subroutine tfm_pre_surfdens_medley2020(wind_north, wind_max, humidity, &
    & accumulation, temperature, surface_density)
    implicit none
    
    real(prec), intent(in) :: wind_north   ! m / s
    real(prec), intent(in) :: wind_max     ! m / s
    real(prec), intent(in) :: humidity     ! ???
    real(prec), intent(in) :: accumulation ! m weq. /s
    real(prec), intent(in) :: temperature  ! K

    real(prec), intent(inout) :: surface_density

    real(prec) :: comp_accumulation

    comp_accumulation = (              &
    &  accumulation                    &
    &  * SECONDS_YEAR                  &
    &  * (WATER_DENSITY / ICE_DENSITY) &
    )

    surface_density = (                &
    &  - (369.6)                       &
    &  + (1.985   * wind_north)        &
    &  + (3.009   * wind_max)          &
    &  - (1.392e5 * humidity)          &
    &  + (27.57   * comp_accumulation) &
    &  + (3.192   * temperature)       &
    )
  end subroutine tfm_pre_surfdens_medley2020


  subroutine tfm_pre_surfdens_ligtenberg2013(surface_temperature, &
    & accumulation_rate, wind_speed_ten_meter, surface_density)
    implicit none

    real(prec), intent(in) :: surface_temperature  ! K
    real(prec), intent(in) :: accumulation_rate    ! m weq. / s
    real(prec), intent(in) :: wind_speed_ten_meter ! m / s

    real(prec), intent(inout) :: surface_density

    real(prec) :: comp_accumulation_rate

    comp_accumulation_rate = (accumulation_rate * 1.0e3 * SECONDS_YEAR)

    surface_density = (                      &
    &  - (151.94)                            &
    &  + (1.4266 * (                         &
    &      (73.6)                            &
    &    + (1.06   * surface_temperature)    &
    &    + (0.0669 * comp_accumulation_rate) &
    &    + (4.77   * wind_speed_ten_meter)   &
    &  ))                                    &
    )
  end subroutine tfm_pre_surfdens_ligtenberg2013
end module tfm_preprocessing
