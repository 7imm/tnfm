module tfm_liquid
  use settings
  use tfm_constants
  implicit none

  real(prec), parameter :: IMP_DENSITY = 830.0

  contains

  subroutine tfm_liquid_bucket(nz, dt, depth, density, temperature, &
    & liquid_water, infiltration_rate, runoff)
    implicit none

    integer, intent(in)                   :: nz
    real(prec), intent(in)                :: dt
    real(prec), dimension(nz), intent(in) :: depth
    real(prec), intent(in)                :: infiltration_rate

    real(prec), dimension(nz), intent(inout) :: density
    real(prec), dimension(nz), intent(inout) :: temperature
    real(prec), dimension(nz), intent(inout) :: liquid_water
    real(prec), intent(inout)               :: runoff

    integer    :: n
    real(prec) :: water
    real(prec) :: dz
    real(prec) :: irr_water_content
    real(prec) :: ice_cap, refreeze_cap, imp_cap
    real(prec) :: storage

    ! infiltrated water mass
    if (infiltration_rate <= 0.0) RETURN
    water = dt * infiltration_rate * WATER_DENSITY

    do n = nz - 1, 1, -1

      dz = (depth(n+1) - depth(n))

      ! potential mass per square meter that might be frozen
      refreeze_cap = (                                                        &
      &  (SPECIFIC_HEAT_ICE * density(n) * dz * (MELT_TEMP - temperature(n))) &
      &  / LATENT_HEAT                                                        &
      )

      ! pore space in kg ice equivalent per square meter
      ice_cap = (1.0 - (density(n) / ICE_DENSITY)) * ICE_DENSITY * dz

      ! pore space available until a impermeable layer is formed
      imp_cap = max(0.0, (IMP_DENSITY - density(n)) * dz)

      ! maximum amount of water to be refrozen
      storage = min(refreeze_cap, ice_cap, water, imp_cap)

      if ( storage < 0.0 ) then
        print *, 'module: tfm_liquid                                   '
        print *, 'subroutine: tfm_liquid_bucket                        '
        print *, 'Variable storage became negative.                    '
        print *, 'strorage: ', storage, 'kg m-2'
        print *, 'Phyiscally this is not possible. Most certainly there'
        print *, 'are temperatures above the melting point occuring.   '
        print *, 'max(temperature): ', maxval(temperature), 'K'
        print *, 'Stopping the simulation right here.                  '
        STOP
      end if

      ! temperature change due to refreezing
      temperature(n) = (                                                   &
      &  temperature(n)                                                    &
      &  + (LATENT_HEAT / (SPECIFIC_HEAT_ICE * density(n) * dz)) * storage &
      )

      ! remaining water
      water = water - storage

      ! density change
      density(n) = density(n) + (storage / dz)
      if ( density(n) >= IMP_DENSITY .or. water <= 0.0 ) EXIT

      ! irreducable water content
      call tfm_liquid_Coleou1998(density(n), irr_water_content)

      ! liquid water remaining in the layer
      storage = min(water, (irr_water_content * WATER_DENSITY * dz))
      liquid_water(n) = liquid_water(n) + storage

      ! remaining water
      water = water - storage
      if ( water <= 0.0 ) EXIT
    end do

    runoff = water ! total runoff water kg
  end subroutine tfm_liquid_bucket


  subroutine tfm_liquid_Coleou1998(density, irr_water_content)
    implicit none

    real(prec), intent(in)    :: density
    real(prec), intent(inout) :: irr_water_content

    irr_water_content = 0.017 + 0.057 * ((ICE_DENSITY - density) / density)
  end subroutine tfm_liquid_Coleou1998
end module tfm_liquid
