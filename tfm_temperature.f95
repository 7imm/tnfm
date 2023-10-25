module tfm_temperature
  use tfm_essentials
  implicit none


  interface
    function dry_thermcond_inter(nz, density, temperature)
      use tfm_essentials
      implicit none

      integer, intent(in)                   :: nz
      real(prec), dimension(nz), intent(in) :: density
      real(prec), dimension(nz), intent(in) :: temperature
      real(prec), dimension(nz)             :: dry_thermcond_inter
    end function dry_thermcond_inter
  end interface


  interface
    function sat_thermcond_inter(nz, density)
      use tfm_essentials
      implicit none

      integer, intent(in)                   :: nz
      real(prec), dimension(nz), intent(in) :: density
      real(prec), dimension(nz)             :: sat_thermcond_inter
    end function sat_thermcond_inter
  end interface


  contains


  function tfm_temperature_diffusion(nz, dt, depth, density, temperature, &
    & heat_capacity, thermal_conductivity) result(d_temperature)
    implicit none

    integer, intent(in)                   :: nz
    real(prec), intent(in)                :: dt
    real(prec), dimension(nz), intent(in) :: depth
    real(prec), dimension(nz), intent(in) :: density
    real(prec), dimension(nz), intent(in) :: temperature
    real(prec), dimension(nz), intent(in) :: heat_capacity
    real(prec), dimension(nz), intent(in) :: thermal_conductivity

    real(prec), dimension(nz) :: d_temperature

    integer                     :: n
    real(prec)                  :: w
    real(prec), dimension(nz-1) :: dz
    real(prec), dimension(nz-2) :: dzt, dzb, dzc
    real(prec), dimension(nz-2) :: kp, kt, kb, tkt, tkb
    real(prec), dimension(nz)   :: q, ap
    real(prec), dimension(nz-1) :: at, ab
    real(prec), dimension(nz)   :: n_temperature

    ! height changes
    dz = depth(2:nz) - depth(1:nz-1)
    dzt = dz(2:nz-1)
    dzb = dz(1:nz-2)
    dzc = (0.5 * (dzt + dzb))

    ! interface conductivity
    kp = thermal_conductivity(2:nz-1)
    tkt = thermal_conductivity(3:nz)
    tkb = thermal_conductivity(1:nz-2)
    kt = (2.0 * kp * tkt) / (kp + tkt)
    kb = (2.0 * kp * tkb) / (kp + tkb)

    ! coefficients
    at(2:nz-1) = -(kt / (dzt * dzc)) * (dt / (density(2:nz-1) * heat_capacity(2:nz-1)))
    ab(1:nz-2) = -(kb / (dzb * dzc)) * (dt / (density(2:nz-1) * heat_capacity(2:nz-1)))
    ap(2:nz-1) = 1.0 - (at(2:nz-1) + ab(1:nz-2))

    ! right hand side
    q = 1.0 * temperature

    ! upper boundary condition
    ap(1) = 1.0
    at(1) = 0.0

    ! lower boundary condition
    ap(nz) = 1.0
    ab(nz-1) = 0.0

    ! Gauß elemination
    do n = 2, nz, 1
      w = ab(n-1) / ap(n-1)
      ap(n) = ap(n) - (w * at(n-1))
      q(n) = q(n) - (w * q(n-1))
    end do

    ! substitution
    n_temperature(nz) = q(nz) / ap(nz)
    do n = nz - 1, 1, -1
      n_temperature(n) = (q(n) - (at(n) * n_temperature(n + 1))) / ap(n)
    end do
    d_temperature = (n_temperature - temperature)
  end function tfm_temperature_diffusion



  function tfm_temperature_airCapacity_CRC(nz, temperature) &
    & result(air_heat_capacity)
    implicit none

    integer, intent(in)                   :: nz
    real(prec), dimension(nz), intent(in) :: temperature
    real(prec), dimension(nz)             :: air_heat_capacity

    integer :: n, m

    real(prec), dimension(25), parameter :: REF_TEMP = (/ &
    &  60.0_prec,                                         &
    &  78.79_prec,                                        &
    &  81.61_prec,                                        &
    &  100.0_prec,                                        &
    &  120.0_prec,                                        &
    &  140.0_prec,                                        &
    &  160.0_prec,                                        &
    &  180.0_prec,                                        &
    &  200.0_prec,                                        &
    &  220.0_prec,                                        &
    &  240.0_prec,                                        &
    &  260.0_prec,                                        &
    &  280.0_prec,                                        &
    &  300.0_prec,                                        &
    &  320.0_prec,                                        &
    &  340.0_prec,                                        &
    &  360.0_prec,                                        &
    &  380.0_prec,                                        &
    &  400.0_prec,                                        &
    &  500.0_prec,                                        &
    &  600.0_prec,                                        &
    &  700.0_prec,                                        &
    &  800.0_prec,                                        &
    &  900.0_prec,                                        &
    &  1000.0_prec                                        &
    /)
    real(prec), dimension(25), parameter :: REF_CAP = (/ &
    &  1901.0_prec,                                      &
    &  1933.0_prec,                                      &
    &  1089.0_prec,                                      &
    &  1040.0_prec,                                      &
    &  1022.0_prec,                                      &
    &  1014.0_prec,                                      &
    &  1011.0_prec,                                      &
    &  1008.0_prec,                                      &
    &  1007.0_prec,                                      &
    &  1006.0_prec,                                      &
    &  1006.0_prec,                                      &
    &  1006.0_prec,                                      &
    &  1006.0_prec,                                      &
    &  1007.0_prec,                                      &
    &  1007.0_prec,                                      &
    &  1009.0_prec,                                      &
    &  1010.0_prec,                                      &
    &  1012.0_prec,                                      &
    &  1014.0_prec,                                      &
    &  1030.0_prec,                                      &
    &  1051.0_prec,                                      &
    &  1075.0_prec,                                      &
    &  1099.0_prec,                                      &
    &  1121.0_prec,                                      &
    &  1141.0_prec                                       &
    /)
! ----------------------------------------------------------------------
! Function: tfm_temperature_airCapacity_CRC
!
! Heat capacity of air at given temperature and pressure of p = 0.1 MPa,
! computed by linear interplation from the vlaues listed in the CRC
! Handbook of Chemistriy and Physics.
!
! Haynes, W. M., Lide, D. R., and Bruno, T. J. (Editors) (2016). Fluid
! Properties. In: CRC Handbook of Chemistry and Physics. CRC Press.
!
! Author: Timm Schultz
!
! Argument:
!   Temperature (K).
!
! Result:
!   air_heat_capacity: Heat capacity of air at given temperature and 
!   p = 0.1 MPa (J kg**-1 K**-1).
! ----------------------------------------------------------------------

    ! loop over all temperature values
    do n = 1, nz, 1

      ! find relevant values for the interpolation
      do m = 1, size(REF_TEMP), 1
        if (REF_TEMP(m) > temperature(n)) then
          EXIT
        end if
      end do
      
      ! linear interpolation
      air_heat_capacity = (                                              &
      &  REF_CAP(m-1)                                                    &
      &  + (                                                             &
      &    ((REF_CAP(m) - REF_CAP(m-1)) / (REF_TEMP(m) - REF_TEMP(m-1))) &
      &    * (temperature(n) - REF_TEMP(m-1))                            &
      &  )                                                               &
      )
    end do
  end function tfm_temperature_airCapacity_CRC


  function tfm_temperature_capacity_paterson1994(nz, density, temperature, &
    & liquid_water) result(n_heat_capacity)
    implicit none
   
    integer, intent(in)                   :: nz
    real(prec), dimension(nz), intent(in) :: density
    real(prec), dimension(nz), intent(in) :: temperature
    real(prec), dimension(nz), intent(in) :: liquid_water
    real(prec), dimension(nz)             :: n_heat_capacity

    real(prec), parameter :: WATER_HEAT_CAPACITY_273 = 4219.4_prec

    real(prec), dimension(nz) :: ice_cap
    real(prec), dimension(nz) :: air_cap
    real(prec), dimension(nz) :: wat_cap
! ----------------------------------------------------------------------
! Function tfm_temperature_capacity_Paterson1994
!
! Heat capacity of firn following the mixture approach described in
! Kaviany (1991) (Voigt model weighting) and using the constant heat 
! capacity of ice given by Paterson (1994).
! The heat capacity of water is that at 273.16 K and is taken from the
! CRC Handbook of Chemistry and Physics (Haynes, 2016).
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of variables "density", "temperature", and
!     "liquid_water".
!   density: Density along the firn profile (kg m**-3).
!   temperature: Temperature along the firn profile (K).
!   liquid_water: Volumetric liquid water content (1).
!
! Result:
!   n_heat_capacity: Effective heat capacity of firn (J kg**-1 K**-1).
! ----------------------------------------------------------------------

    ice_cap = 2009.0_prec
    wat_cap = WATER_HEAT_CAPACITY_273
    air_cap = tfm_temperature_airCapacity_CRC(nz, temperature)

    n_heat_capacity = (                                                  &
    &  (ice_cap * (density / ICE_DENSITY))                               &
    &  + (wat_cap * liquid_water)                                        &
    &  + (air_cap * (1.0_prec - (density / ICE_DENSITY) - liquid_water)) &
    )
  end function tfm_temperature_capacity_paterson1994


  function tfm_temperature_capacity_Cuffey2010(nz, density, temperature, &
    & liquid_water) result(n_heat_capacity)
    implicit none

    integer, intent(in)                   :: nz
    real(prec), dimension(nz), intent(in) :: density
    real(prec), dimension(nz), intent(in) :: temperature
    real(prec), dimension(nz), intent(in) :: liquid_water
    real(prec), dimension(nz)             :: n_heat_capacity

    real(prec), parameter :: WATER_HEAT_CAPACITY_273 = 4219.4_prec

    real(prec), dimension(nz) :: ice_cap
    real(prec), dimension(nz) :: wat_cap
    real(prec), dimension(nz) :: air_cap
! ----------------------------------------------------------------------
! Function tfm_temperature_capacity_Cuffey2010
!
! Heat capacity of firn following the mixture approach described in
! Kaviany (1991) (Voigt model weighting) and using the temperature
! dependentheat capacity of ice given by Cuffey & Paterson (2010).
! The heat capacity of water is that at 273.16 K and is taken from the
! CRC Handbook of Chemistry and Physics (Haynes, 2016).
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of variables "density", "temperature", and
!     "liquid_water".
!   density: Density along the firn profile (kg m**-3).
!   temperature: Temperature along the firn profile (K).
!   liquid_water: Volumetric liquid water content (1).
!
! Result:
!   n_heat_capacity: Effective heat capacity of firn (J kg**-1 K**-1).
! ----------------------------------------------------------------------

    ice_cap = (152.5_prec + (7.122_prec * temperature))
    wat_cap = WATER_HEAT_CAPACITY_273
    air_cap = tfm_temperature_airCapacity_CRC(nz, temperature)

    n_heat_capacity = (                                                  &
    &  (ice_cap * (density / ICE_DENSITY))                               &
    &  + (wat_cap * liquid_water)                                        &
    &  + (air_cap * (1.0_prec - (density / ICE_DENSITY) - liquid_water)) &
    )
  end function tfm_temperature_capacity_Cuffey2010


  function tfm_temperature_airConductivity_CRC(nz, temperature) &
    & result(air_thermal_conductivity)
    implicit none

    integer, intent(in)                   :: nz
    real(prec), dimension(nz), intent(in) :: temperature
    real(prec), dimension(nz)             :: air_thermal_conductivity

    integer :: n, m

    real(prec), dimension(25), parameter :: REF_TEMP = (/ &
    &  60.0_prec,                                         &
    &  78.79_prec,                                        &
    &  81.61_prec,                                        &
    &  100.0_prec,                                        &
    &  120.0_prec,                                        &
    &  140.0_prec,                                        &
    &  160.0_prec,                                        &
    &  180.0_prec,                                        &
    &  200.0_prec,                                        &
    &  220.0_prec,                                        &
    &  240.0_prec,                                        &
    &  260.0_prec,                                        &
    &  280.0_prec,                                        &
    &  300.0_prec,                                        &
    &  320.0_prec,                                        &
    &  340.0_prec,                                        &
    &  360.0_prec,                                        &
    &  380.0_prec,                                        &
    &  400.0_prec,                                        &
    &  500.0_prec,                                        &
    &  600.0_prec,                                        &
    &  700.0_prec,                                        &
    &  800.0_prec,                                        &
    &  900.0_prec,                                        &
    &  1000.0_prec                                        &
    /)
    real(prec), dimension(25), parameter :: REF_COND = (/ &
    &  0.1711_prec,                                       &
    &  0.1402_prec,                                       &
    &  0.7673_prec,                                       &
    &  0.9469_prec,                                       &
    &  1.138_prec,                                        &
    &  1.324_prec,                                        &
    &  1.505_prec,                                        &
    &  1.68_prec,                                         &
    &  1.85_prec,                                         &
    &  2.016_prec,                                        &
    &  2.177_prec,                                        &
    &  2.335_prec,                                        &
    &  2.488_prec,                                        &
    &  2.638_prec,                                        &
    &  2.785_prec,                                        &
    &  2.929_prec,                                        &
    &  3.071_prec,                                        &
    &  3.209_prec,                                        &
    &  3.345_prec,                                        &
    &  3.994_prec,                                        &
    &  4.601_prec,                                        &
    &  5.176_prec,                                        &
    &  5.725_prec,                                        &
    &  6.254_prec,                                        &
    &  6.768_prec                                         &
    /)
! ----------------------------------------------------------------------
! Function: tfm_temperature_airConductivity_CRC
!
! Thermal conductivity at given temperature and pressure of p = 0.1 MPa,
! computed by linear interpolation from the values listed in the CRC
! Handbook of Chemistry and Physics.
!
! Haynes, W. M., Lide, D. R., and Bruno, T. J. (Editors) (2016). Fluid
! Properties. In: CRC Handbook of Chemistry and Physics. CRC Press.
!
! Author: Timm Schultz
!
! Argument:
!   Temperature (K).
!
! Result:
!   air_thermal_conductivity: Thermal conductivity of air at given
!     temperature and p = 0.1 MPa (W m**-1 K**-1).
! ----------------------------------------------------------------------

    ! loop over all temperature values
    do n = 1, nz, 1

      ! find relevant values for the interpolation
      do m = 1, size(REF_TEMP), 1
        if (REF_TEMP(m) > temperature(n)) then
          EXIT
        end if
      end do
      
      ! linear interpolation
      air_thermal_conductivity(n) = (                                      &
      &  REF_COND(m-1)                                                     &
      &  + (                                                               &
      &    ((REF_COND(m) - REF_COND(m-1)) / (REF_TEMP(m) - REF_TEMP(m-1))) &
      &    * (temperature(n) - REF_TEMP(m-1))                              &
      &  )                                                                 &
      )
    end do
  end function tfm_temperature_airConductivity_CRC


  function tfm_temperature_iceConductivity_Cuffey2010(nz, temperature) &
    & result(thermal_conductivity)
    implicit none

    integer, intent(in)                   :: nz
    real(prec), dimension(nz), intent(in) :: temperature
    real(prec), dimension(nz)             :: thermal_conductivity
! ----------------------------------------------------------------------
! Function: tfm_temperature_iceConductivity_Cuffey2010
!
! Cuffey, K. M. and Paterson, W. S. B. (2010). The Physics of Glaciers.
! p. 400.
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimenstion of variable "temperature".
!   temperature: Temperature (K).
!
! Result:
!   thermal_conductivity: Thermal conductivity of ice at given
!     temperature (W m**-1 K**-1).
! ----------------------------------------------------------------------

   thermal_conductivity = (9.828_prec * exp(-4.7E-3_prec * temperature))
  end function tfm_temperature_iceConductivity_Cuffey2010


  function tfm_temperature_mixture_miller1969UpperBound(nz, porosity, &
    & solid_cond, liquid_cond) result(eff_cond)
    implicit none

    integer, intent(in)                   :: nz
    real(prec), dimension(nz), intent(in) :: porosity 
    real(prec), dimension(nz), intent(in) :: solid_cond
    real(prec), dimension(nz), intent(in) :: liquid_cond
    real(prec), dimension(nz)             :: eff_cond

    real(prec), parameter     :: PARAM_G = (1.0_prec / 9.0_prec)
    real(prec), dimension(nz) :: fu
! ----------------------------------------------------------------------
! Function: tfm_temperature_mixutre_miller1969UpperBound
!
! Miller, M. N. (1969). Bounds for Effective Electrical, Thermal, and
! Magenetic Properties of Heterogeneous Materials. Journal of
! Mathematical Physics, Volume 10, Number 11, pp. 1988-2004,
! https://doi.org/10.1063/1.1664794
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of variabels "density" and "temperature".
!   porosity: Porosity (1).
!   solid_cond: Thermal conductivity of the solid phase (W m**-1 K**-1).
!   liquid_cond: Thermal conductvitiy of the liquid phase
!     (W m**-1 K**-1).
!
! Result:
!   eff_cond: Effective thermal conductivity (W m**-1 K**-1).
! ----------------------------------------------------------------------
    fu = (1.0_prec - (                                                          &
    &  (                                                                        &
    &    (1.0_prec - porosity)                                                  &
    &    * (((solid_cond / liquid_cond) - 1.0_prec)**2.0_prec)                  &
    &    * porosity                                                             &
    &  )                                                                        &
    &  / (                                                                      &
    &    3.0_prec                                                               &
    &    * (                                                                    &
    &      1.0_prec                                                             &
    &      + ((1.0_prec - porosity) * ((solid_cond / liquid_cond) - 1.0_prec))  &
    &    )                                                                      &
    &    * (                                                                    &
    &      1.0_prec                                                             &
    &      + (                                                                  &
    &        (1.0_prec - porosity)                                              &
    &        + ((3.0_prec * ((2.0_prec * porosity) - 1.0_prec)) * PARAM_G)      &
    &      )                                                                    &
    &      * ((solid_cond / liquid_cond) - 1.0_prec)                            &
    &    )                                                                      &
    &  )                                                                        &
    ))

    eff_cond = (fu * (                                                     &
    &  1.0_prec                                                            &
    &  + ((1.0_prec - porosity) * ((solid_cond / liquid_cond) - 1.0_prec)) &
    ))
  end function tfm_temperature_mixture_miller1969UpperBound


  function tfm_temperature_mixture_miller1969LowerBound(nz, porosity, &
    & solid_cond, liquid_cond) result(eff_cond)
    implicit none

    integer, intent(in)                   :: nz
    real(prec), dimension(nz), intent(in) :: porosity
    real(prec), dimension(nz), intent(in) :: solid_cond
    real(prec), dimension(nz), intent(in) :: liquid_cond
    real(prec), dimension(nz)             :: eff_cond

    real(prec), parameter     :: PARAM_G = (1.0_prec / 9.0_prec)
! ----------------------------------------------------------------------
! Function: tfm_temperature_mixutre_miller1969LowerBound
!
! Miller, M. N. (1969). Bounds for Effective Electrical, Thermal, and
! Magenetic Properties of Heterogeneous Materials. Journal of
! Mathematical Physics, Volume 10, Number 11, pp. 1988-2004,
! https://doi.org/10.1063/1.1664794
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of variabels "density" and "temperature".
!   porosity: Porosity (1).
!   solid_cond: Thermal conductivity of the solid phase (W m**-1 K**-1).
!   liquid_cond: Thermal conductvitiy of the liquid phase
!     (W m**-1 K**-1).
!
! Result:
!   eff_cond: Effective thermal conductivity (W m**-1 K**-1).
! ----------------------------------------------------------------------

    eff_cond = (                                                             &
    &  solid_cond * ((                                                       &
    &    (solid_cond / liquid_cond)                                          &
    &    - ((1.0_prec - porosity) * ((solid_cond / liquid_cond) - 1.0_prec)) &
    &    - (                                                                 &
    &      (                                                                 &
    &        (4.0_prec / 3.0_prec)                                           &
    &        * (((solid_cond / liquid_cond) - 1.0_prec)**2.0_prec)           &
    &        * (1.0_prec - porosity)                                         &
    &        * porosity                                                      &
    &      )                                                                 &
    &      / (                                                               &
    &        1.0_prec                                                        &
    &        + (solid_cond / liquid_cond)                                    &
    &        + (                                                             &
    &          (3.0_prec * (1.0_prec - (2.0_prec * porosity)))               &
    &          * ((solid_cond / liquid_cond) - 1.0_prec)                     &
    &          * PARAM_G                                                     &
    &        )                                                               &
    &      )                                                                 &
    &    )                                                                   &
    &  )**(-1.0_prec))                                                       &
    )
  end function tfm_temperature_mixture_miller1969LowerBound


  function tfm_temperature_mixture_geomMean(nz, porosity, solid_cond, &
    & liquid_cond) result(eff_cond)
    implicit none

    integer, intent(in)                   :: nz
    real(prec), dimension(nz), intent(in) :: porosity
    real(prec), dimension(nz), intent(in) :: solid_cond
    real(prec), dimension(nz), intent(in) :: liquid_cond
    real(prec), dimension(nz)             :: eff_cond
! ----------------------------------------------------------------------
! Function: tfm_temperature_mixture_geomMean
!
! See:
! Kaviany, M. (1991). Principles of Heat Transfer in Porous Media. in:
! Ling, F. F. (edt.) Mechanical Engineering Series, Springer-Verlag, New
! York, Berlin, Heidelberg, https://doi.org/10.1063/1.1664794
! p. 126
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of variabels "density" and "temperature".
!   porosity: Porosity (1).
!   solid_cond: Thermal conductivity of the solid phase (W m**-1 K**-1).
!   liquid_cond: Thermal conductvitiy of the liquid phase
!     (W m**-1 K**-1).
!
! Result:
!   eff_cond: Effective thermal conductivity (W m**-1 K**-1).
! ----------------------------------------------------------------------

    eff_cond = (                              &
    &  (liquid_cond**porosity)                &
    &  * (solid_cond**(1.0_prec - porosity))  &
    )
  end function tfm_temperature_mixture_geomMean


  function tfm_temperature_mixture_voigt(nz, porosity, solid_cond, &
    & liquid_cond) result(eff_cond)
    implicit none

    integer, intent(in)                   :: nz
    real(prec), dimension(nz), intent(in) :: porosity
    real(prec), dimension(nz), intent(in) :: solid_cond
    real(prec), dimension(nz), intent(in) :: liquid_cond
    real(prec), dimension(nz)             :: eff_cond
! ----------------------------------------------------------------------
! Function tfm_temperature_mixture_voigt
!
! Voigt, W. (1889). Ueber die Beziehung zwischen den beiden
! Elasticitätsconstanten Isotroper Körper. Annalen der Physik, Volume
! 274, Issue 12, pp. 573-587,
! https://doi.org/10.1002/andp.18892741206
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of variabels "density" and "temperature".
!   porosity: Porosity (1).
!   solid_cond: Thermal conductivity of the solid phase (W m**-1 K**-1).
!   liquid_cond: Thermal conductvitiy of the liquid phase
!     (W m**-1 K**-1).
!
! Result:
!   eff_cond: Effective thermal conductivity (W m**-1 K**-1).
! ----------------------------------------------------------------------
    
    eff_cond = (                              &
    &  (porosity * liquid_cond)               &
    &  + ((1.0_prec - porosity) * solid_cond) &
    )
  end function tfm_temperature_mixture_voigt


  function tfm_temperature_mixture_reuss(nz, porosity, solid_cond, &
    & liquid_cond) result(eff_cond)
    implicit none

    integer, intent(in)                   :: nz
    real(prec), dimension(nz), intent(in) :: porosity
    real(prec), dimension(nz), intent(in) :: solid_cond
    real(prec), dimension(nz), intent(in) :: liquid_cond
    real(prec), dimension(nz)             :: eff_cond
! ----------------------------------------------------------------------
! Function tfm_temperature_mixture_reuss
!
! Reuss, A. (1929). Berechnung der Fließgrenze von Mischkristallen auf
! Grund der Platizitätsbedingung für Einkristalle. ZAMM - Journal of
! Applied Mathematics and Mechanics, Volume 9, Issue 1, pp. 49-58,
! https://doi.org/10.1002/zamm.19290090104
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of variabels "density" and "temperature".
!   porosity: Porosity (1).
!   solid_cond: Thermal conductivity of the solid phase (W m**-1 K**-1).
!   liquid_cond: Thermal conductvitiy of the liquid phase
!     (W m**-1 K**-1).
!
! Result:
!   eff_cond: Effective thermal conductivity (W m**-1 K**-1).
! ----------------------------------------------------------------------
    
    eff_cond = (                                                         &
    &  solid_cond                                                        &
    &  / ((porosity * (solid_cond / liquid_cond)) + 1.0_prec - porosity) &
    )
  end function tfm_temperature_mixture_reuss


  function tfm_temperature_conduct_calonne2019(nz, density, temperature) &
    & result(n_thermal_conductivity)
    implicit none

    integer, intent(in) :: nz
    real(prec), dimension(nz), intent(in) :: density
    real(prec), dimension(nz), intent(in) :: temperature
    real(prec), dimension(nz)             :: n_thermal_conductivity

    real(prec), dimension(nz) :: theta
    real(prec), dimension(nz) :: firn_cond
    real(prec), dimension(nz) :: snow_cond
    real(prec), dimension(nz) :: air_cond
    real(prec), dimension(nz) :: ice_cond

    real(prec), parameter :: TRANSITION_DENSITY = 450.0_prec ! kg m**-3
    real(prec), parameter :: ICE_COND_REF = 2.107_prec ! W m**-1 K**-1
    real(prec), parameter :: AIR_COND_REF = 0.024_prec ! W m**-1 K**-1
    real(prec), parameter :: PARAM_A = 0.02 ! m**3 kf**-1
! ----------------------------------------------------------------------
! Function: tfm_temperature_conduct_calonne2019
!
! Calonne, N., Milliancourt, L., Burr, A., Philip, A., Martin, C. L.,
! Flin, F., and Geindreau, C. (2019). Thermal Conductivity of Snow,
! Firn, and Porous Ice From 3-D Image-Based Computations, Volume 46,
! Issue 22, pp. 13079-13089, https://doi.org/10.1029/2019GL085228
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of variabels "density" and "temperature".
!   density: Density along the firn profile (kg m**-3).
!   temperature: Temperature along the firn profile (K).
!
! Result:
!   n_thermal_conductivity: Thermal conductvitiy along the firn profile
!     (W m**-1 K**-1).
! ----------------------------------------------------------------------
    
    theta = (                                                        &
    &  1.0_prec                                                      &
    &  / (                                                           &
    &    1.0_prec                                                    &
    &    + exp(-2.0_prec * PARAM_A * (density - TRANSITION_DENSITY)) &
    &  )                                                             &
    )

    snow_cond = (                            &
    &  0.024_prec                            &
    &  - (1.23E-4_prec * density)            &
    &  + (2.4E-6_prec * (density**2.0_prec)) &
    )

    firn_cond = (                                  &
    &  2.107_prec                                  &
    &  + (0.003618_prec * (density - ICE_DENSITY)) &
    )

    air_cond = tfm_temperature_airConductivity_CRC(nz, temperature)
    ice_cond = tfm_temperature_iceConductivity_Cuffey2010(nz, temperature)

    n_thermal_conductivity = (                                   &
    &  (                                                         &
    &    (1.0_prec - theta) * snow_cond * (                      &
    &      (ice_cond * air_cond) / (ICE_COND_REF * AIR_COND_REF) &
    &    )                                                       &
    &  )                                                         &
    &  + (theta * firn_cond * (ice_cond / ICE_COND_REF))         &
    )
  end function tfm_temperature_conduct_calonne2019


  function tfm_temperature_conduct_marchenko2019(nz, density, temperature) &
    & result(n_thermal_conductivity)
    implicit none

    integer, intent(in)                   :: nz
    real(prec), dimension(nz), intent(in) :: density
    real(prec), dimension(nz), intent(in) :: temperature
    real(prec), dimension(nz)             :: n_thermal_conductivity
! ----------------------------------------------------------------------
! Function: tfm_temperature_conduct_marchenko2019
!
! Marchenko, S., Cheng, G., Lötstedt, P., Pohjola, V., Pettersson, R.,
! van Pelt, W., and Reijmer, C. (2019). Thermal conductivity of firn at
! Lomonosovfonna, Svalbard, derived from subsurface temperature
! measurements. The Cryosphere, Volume 13, Issue 7, pp. 1843-1859,
! https://doi.org/10.5194/tc-13-1843-2019
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of variabels "density" and "temperature".
!   density: Density along the firn profile (kg m**-3).
!   temperature: Temperature along the firn profile (K).
!
! Result:
!   n_thermal_conductivity: Thermal conductvitiy along the firn profile
!     (W m**-1 K**-1).
! ----------------------------------------------------------------------
    
    call tfm_essentials_do_nothing(nz, temperature)
    
    n_thermal_conductivity = ((0301E-2_prec * density) - 0.724_prec)
  end function tfm_temperature_conduct_marchenko2019


  function tfm_temperature_conduct_sturm1997(nz, density, temperature) &
    & result(n_thermal_conductivity)
    implicit none

    integer, intent(in)                   :: nz
    real(prec), dimension(nz), intent(in) :: density
    real(prec), dimension(nz), intent(in) :: temperature
    real(prec), dimension(nz)             :: n_thermal_conductivity
! ----------------------------------------------------------------------
! Function: tfm_temperature_conduct_sturm1997
!
! Sturm, M., Holmgren, J., König, M., and Morris, K. (1997). The thermal
! conductivity of seasonal snow. Journal of Glaciology, Volume 43, Issue
! 143, pp. 26-41, https://doi.org/10.3189/S0022143000002781
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of variabels "density" and "temperature".
!   density: Density along the firn profile (kg m**-3).
!   temperature: Temperature along the firn profile (K).
!
! Result:
!   n_thermal_conductivity: Thermal conductvitiy along the firn profile
!     (W m**-1 K**-1).
! ----------------------------------------------------------------------
    
    call tfm_essentials_do_nothing(nz, temperature)

    n_thermal_conductivity = (         &
    &  (0.138)                         &
    &  - ((1.010e-3) * density)        &
    &  + ((3.233e-6) * (density**2.0)) &
    )
  end function tfm_temperature_conduct_sturm1997


  function tfm_temperature_conduct_miller1969UpperBound(nz, density, temperature) &
    & result(n_thermal_conductivity)
    implicit none

    integer, intent(in)                   :: nz
    real(prec), dimension(nz), intent(in) :: density
    real(prec), dimension(nz), intent(in) :: temperature
    real(prec), dimension(nz)             :: n_thermal_conductivity

    real(prec), dimension(nz) :: air_cond
    real(prec), dimension(nz) :: ice_cond
    real(prec), dimension(nz) :: porosity
! ----------------------------------------------------------------------
! Function: tfm_temperature_conduct_miller1969UpperBound
!
! Miller, M. N. (1969). Bounds for Effective Electrical, Thermal, and
! Magenetic Properties of Heterogeneous Materials. Journal of
! Mathematical Physics, Volume 10, Number 11, pp. 1988-2004,
! https://doi.org/10.1063/1.1664794
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of variabels "density" and "temperature".
!   density: Density along the firn profile (kg m**-3).
!   temperature: Temperature along the firn profile (K).
!
! Result:
!   n_thermal_conductivity: Thermal conductvitiy along the firn profile
!     (W m**-1 K**-1).
! ----------------------------------------------------------------------
    
    air_cond = tfm_temperature_airConductivity_CRC(nz, temperature)
    ice_cond = tfm_temperature_iceConductivity_Cuffey2010(nz, temperature)
    porosity = (1.0_prec - (density / ICE_DENSITY))

    n_thermal_conductivity = tfm_temperature_mixture_miller1969UpperBound( &
      nz, porosity, ice_cond, air_cond                                     &
    )
  end function tfm_temperature_conduct_miller1969UpperBound


  function tfm_temperature_conduct_miller1969LowerBound(nz, density, temperature) &
    & result(n_thermal_conductivity)
    implicit none

    integer, intent(in)                   :: nz
    real(prec), dimension(nz), intent(in) :: density
    real(prec), dimension(nz), intent(in) :: temperature
    real(prec), dimension(nz)             :: n_thermal_conductivity

    real(prec), parameter     :: PARAM_G = (1.0_prec / 9.0_prec)
    real(prec), dimension(nz) :: air_cond
    real(prec), dimension(nz) :: ice_cond
    real(prec), dimension(nz) :: porosity
! ----------------------------------------------------------------------
! Function: tfm_temperature_conduct_miller1969LowerBound
!
! Miller, M. N. (1969). Bounds for Effective Electrical, Thermal, and
! Magenetic Properties of Heterogeneous Materials. Journal of
! Mathematical Physics, Volume 10, Number 11, pp. 1988-2004,
! https://doi.org/10.1063/1.1664794
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of variabels "density" and "temperature".
!   porosity: Porosity (1).
!   solid_cond: Thermal conductivity of the solid phase (W m**-1 K**-1).
!   liquid_cond: Thermal conductvitiy of the liquid phase
!     (W m**-1 K**-1).
!
! Result:
!   eff_cond: Effective thermal conductivity (W m**-1 K**-1).
! ----------------------------------------------------------------------

    air_cond = tfm_temperature_airConductivity_CRC(nz, temperature)
    ice_cond = tfm_temperature_iceConductivity_Cuffey2010(nz, temperature)
    porosity = (1.0_prec - (density / ICE_DENSITY))

    n_thermal_conductivity = tfm_temperature_mixture_miller1969LowerBound( &
      nz, porosity, ice_cond, air_cond                                     &
    )
  end function tfm_temperature_conduct_miller1969LowerBound

  
  function tfm_temperature_conduct_geomMean(nz, density, temperature) &
    & result(n_thermal_conductivity)
    implicit none

    integer, intent(in)                   :: nz
    real(prec), dimension(nz), intent(in) :: density
    real(prec), dimension(nz), intent(in) :: temperature
    real(prec), dimension(nz)             :: n_thermal_conductivity

    real(prec), dimension(nz) :: air_cond
    real(prec), dimension(nz) :: ice_cond
    real(prec), dimension(nz) :: porosity
! ----------------------------------------------------------------------
! Function: tfm_temperature_conduct_geomMean
!
! See:
! Kaviany, M. (1991). Principles of Heat Transfer in Porous Media. in:
! Ling, F. F. (edt.) Mechanical Engineering Series, Springer-Verlag, New
! York, Berlin, Heidelberg, https://doi.org/10.1063/1.1664794
! p. 126
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of variabels "density" and "temperature".
!   density: Density along the firn profile (kg m**-3).
!   temperature: Temperature along the firn profile (K).
!
! Result:
!   n_thermal_conductivity: Thermal conductvitiy along the firn profile
!     (W m**-1 K**-1).
! ----------------------------------------------------------------------

    air_cond = tfm_temperature_airConductivity_CRC(nz, temperature)
    ice_cond = tfm_temperature_iceConductivity_Cuffey2010(nz, temperature)
    porosity = (1.0_prec - (density / ICE_DENSITY))

    n_thermal_conductivity = tfm_temperature_mixture_geomMean( &
    &  nz, porosity, ice_cond, air_cond                        &
    )
  end function tfm_temperature_conduct_geomMean


  function tfm_temperature_sat_cond_Voigt(nz, density) &
    & result (sat_cond)
    implicit none
    
    integer, intent(in)                   :: nz
    real(prec), dimension(nz), intent(in) :: density
    real(prec), dimension(nz)             :: sat_cond

    ! Thermal conductivity of water at 273.15 K (W m**-1 K**-1)
    real(prec), parameter     :: WATER_COND_273 = 0.561_prec
    real(prec), dimension(nz) :: porosity
    real(prec), dimension(nz) :: temperature
    real(prec), dimension(nz) :: ice_cond
    real(prec), dimension(nz) :: water_cond
! ----------------------------------------------------------------------
! Function: tfm_temperature_sat_Voigt
!
! Function calculating the effective thermal conductivity of firn at
! water saturation (needed for the computation of the effective thermal
! conductivity at unsaturated water levels). The effective thermal
! conductivity is computed using a the model of Voigt (1889). The firn
! temperature is assumed to be the melt temperature.
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of variables "density" and "temperature".
!   density: Density along the firn profile (kg m**-3).
!   temperature: Temperature along the firn profile (K).
!
! Result:
!   sat_cond: Effective thermal conductivity of firn at water
!     saturation (W m**-1 K**-1).
! ----------------------------------------------------------------------

    porosity = (1.0_prec - (density / ICE_DENSITY))
    temperature(:) = MELT_TEMP
    ice_cond = tfm_temperature_iceConductivity_Cuffey2010(nz, temperature)
    water_cond(:) = WATER_COND_273

    sat_cond = tfm_temperature_mixture_Voigt(nz, porosity, ice_cond, WATER_COND)
  end function tfm_temperature_sat_cond_Voigt


  function tfm_temperature_sat_cond_Reuss(nz, density) &
    & result (sat_cond)
    implicit none
    
    integer, intent(in)                   :: nz
    real(prec), dimension(nz), intent(in) :: density
    real(prec), dimension(nz)             :: sat_cond

    ! Thermal conductivity of water at 273.15 K (W m**-1 K**-1)
    real(prec), parameter     :: WATER_COND_273 = 0.561_prec
    real(prec), dimension(nz) :: porosity
    real(prec), dimension(nz) :: temperature
    real(prec), dimension(nz) :: ice_cond
    real(prec), dimension(nz) :: water_cond
! ----------------------------------------------------------------------
! Function: tfm_temperature_sat_Reuss
!
! Function calculating the effective thermal conductivity of firn at
! water saturation (needed for the computation of the effective thermal
! conductivity at unsaturated water levels). The effective thermal
! conductivity is computed using a the model of Reuss (1929). The firn
! temperature is assumed to be the melt temperature.
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of variables "density" and "temperature".
!   density: Density along the firn profile (kg m**-3).
!   temperature: Temperature along the firn profile (K).
!
! Result:
!   sat_cond: Effective thermal conductivity of firn at water
!     saturation (W m**-1 K**-1).
! ----------------------------------------------------------------------

    porosity = (1.0_prec - (density / ICE_DENSITY))
    temperature(:) = MELT_TEMP
    ice_cond = tfm_temperature_iceConductivity_Cuffey2010(nz, temperature)
    water_cond(:) = WATER_COND_273

    sat_cond = tfm_temperature_mixture_Reuss(nz, porosity, ice_cond, WATER_COND)
  end function tfm_temperature_sat_cond_Reuss


  function tfm_temperature_sat_cond_Miller1969UpperBound(nz, density) &
    & result (sat_cond)
    implicit none
    
    integer, intent(in)                   :: nz
    real(prec), dimension(nz), intent(in) :: density
    real(prec), dimension(nz)             :: sat_cond

    ! Thermal conductivity of water at 273.15 K (W m**-1 K**-1)
    real(prec), parameter     :: WATER_COND_273 = 0.561_prec
    real(prec), dimension(nz) :: porosity
    real(prec), dimension(nz) :: temperature
    real(prec), dimension(nz) :: ice_cond
    real(prec), dimension(nz) :: water_cond
! ----------------------------------------------------------------------
! Function: tfm_temperature_sat_Miller1969UpperBound
!
! Function calculating the effective thermal conductivity of firn at
! water saturation (needed for the computation of the effective thermal
! conductivity at unsaturated water levels). The effective thermal
! conductivity is computed using a the upper bound model of Miller
! (1969). The firn temperature is assumed to be the melt temperature.
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of variables "density" and "temperature".
!   density: Density along the firn profile (kg m**-3).
!   temperature: Temperature along the firn profile (K).
!
! Result:
!   sat_cond: Effective thermal conductivity of firn at water
!     saturation (W m**-1 K**-1).
! ----------------------------------------------------------------------

    porosity = (1.0_prec - (density / ICE_DENSITY))
    temperature(:) = MELT_TEMP
    ice_cond = tfm_temperature_iceConductivity_Cuffey2010(nz, temperature)
    water_cond(:) = WATER_COND_273

    sat_cond = tfm_temperature_mixture_Miller1969UpperBound( &
    & nz, porosity, ice_cond, WATER_COND                     &
    )
  end function tfm_temperature_sat_cond_Miller1969UpperBound


  function tfm_temperature_sat_cond_Miller1969LowerBound(nz, density) &
    & result (sat_cond)
    implicit none
    
    integer, intent(in)                   :: nz
    real(prec), dimension(nz), intent(in) :: density
    real(prec), dimension(nz)             :: sat_cond

    ! Thermal conductivity of water at 273.15 K (W m**-1 K**-1)
    real(prec), parameter     :: WATER_COND_273 = 0.561_prec
    real(prec), dimension(nz) :: porosity
    real(prec), dimension(nz) :: temperature
    real(prec), dimension(nz) :: ice_cond
    real(prec), dimension(nz) :: water_cond
! ----------------------------------------------------------------------
! Function: tfm_temperature_sat_Miller1969LowerBound
!
! Function calculating the effective thermal conductivity of firn at
! water saturation (needed for the computation of the effective thermal
! conductivity at unsaturated water levels). The effective thermal
! conductivity is computed using a the lower bound model of Miller
! (1969). The firn temperature is assumed to be the melt temperature.
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of variables "density" and "temperature".
!   density: Density along the firn profile (kg m**-3).
!   temperature: Temperature along the firn profile (K).
!
! Result:
!   sat_cond: Effective thermal conductivity of firn at water
!     saturation (W m**-1 K**-1).
! ----------------------------------------------------------------------

    porosity = (1.0_prec - (density / ICE_DENSITY))
    temperature(:) = MELT_TEMP
    ice_cond = tfm_temperature_iceConductivity_Cuffey2010(nz, temperature)
    water_cond(:) = WATER_COND_273

    sat_cond = tfm_temperature_mixture_Miller1969LowerBound( &
    & nz, porosity, ice_cond, WATER_COND                     &
    )
  end function tfm_temperature_sat_cond_Miller1969LowerBound


  function tfm_temperature_sat_cond_geomMean(nz, density) &
    & result (sat_cond)
    implicit none
    
    integer, intent(in)                   :: nz
    real(prec), dimension(nz), intent(in) :: density
    real(prec), dimension(nz)             :: sat_cond

    ! Thermal conductivity of water at 273.15 K (W m**-1 K**-1)
    real(prec), parameter     :: WATER_COND_273 = 0.561_prec
    real(prec), dimension(nz) :: porosity
    real(prec), dimension(nz) :: temperature
    real(prec), dimension(nz) :: ice_cond
    real(prec), dimension(nz) :: water_cond
! ----------------------------------------------------------------------
! Function: tfm_temperature_sat_geomMean
!
! Function calculating the effective thermal conductivity of firn at
! water saturation (needed for the computation of the effective thermal
! conductivity at unsaturated water levels). The effective thermal
! conductivity is computed using a geometric mean weighting approach
! (following Kaviany, 1991). The firn temperature is assumed to be the
! melt temperature.
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of variables "density" and "temperature".
!   density: Density along the firn profile (kg m**-3).
!   temperature: Temperature along the firn profile (K).
!
! Result:
!   sat_cond: Effective thermal conductivity of firn at water
!     saturation (W m**-1 K**-1).
! ----------------------------------------------------------------------

    porosity = (1.0_prec - (density / ICE_DENSITY))
    temperature(:) = MELT_TEMP
    ice_cond = tfm_temperature_iceConductivity_Cuffey2010(nz, temperature)
    water_cond(:) = WATER_COND_273

    sat_cond = tfm_temperature_mixture_geomMean(nz, porosity, ice_cond, WATER_COND)
  end function tfm_temperature_sat_cond_geomMean


  function tfm_temperature_liquid_cond_geomMean(nz, density, temperature, &
    & liquid_water, dry_thermcond_model, sat_thermcond_model)             &
    & result(n_thermal_conductivity)
    implicit none

    integer, intent(in)                     :: nz
    real(prec), dimension(nz), intent(in)   :: density
    real(prec), dimension(nz), intent(in)   :: temperature
    real(prec), dimension(nz), intent(in)   :: liquid_water
    procedure(dry_thermcond_inter), pointer :: dry_thermcond_model
    procedure(sat_thermcond_inter), pointer :: sat_thermcond_model
    real(prec), dimension(nz)               :: n_thermal_conductivity

    real(prec), dimension(nz) :: dry_cond
    real(prec), dimension(nz) :: sat_cond
! ----------------------------------------------------------------------
! Function: tfm_temperature_liquid_cond_geomMean
!
! The function computes the effective thermal conductivity of
! unsaturated firn following Kaviany (1991), using geometric mean
! weighting between dry and saturated conditions.
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of variables "density", "temperature", and
!     "liuqid_water".
!   density: Density along the firn profile (kg m**-3).
!   temperature: Temperature along the firn profile (K).
!   liquid_water: Liquid water content along the profile (1).
!   dry_thermcond_model: Model to use for the computation of effective
!     thermal conductvity at dry conditions.
!   sat_thermcond_model: Model to use for the computation of effective
!     thermal conductviity at saturated conditions.
! ----------------------------------------------------------------------

    dry_cond = dry_thermcond_model(nz, density, temperature)
    sat_cond = sat_thermcond_model(nz, density)

    n_thermal_conductivity = tfm_temperature_mixture_geomMean( &
    & nz, liquid_water, dry_cond, sat_cond                     &
    )
  end function tfm_temperature_liquid_cond_geomMean


  function tfm_temperature_liquid_cond_Voigt(nz, density, temperature, &
    & liquid_water, dry_thermcond_model, sat_thermcond_model)          &
    & result(n_thermal_conductivity)
    implicit none

    integer, intent(in)                     :: nz
    real(prec), dimension(nz), intent(in)   :: density
    real(prec), dimension(nz), intent(in)   :: temperature
    real(prec), dimension(nz), intent(in)   :: liquid_water
    procedure(dry_thermcond_inter), pointer :: dry_thermcond_model
    procedure(sat_thermcond_inter), pointer :: sat_thermcond_model
    real(prec), dimension(nz)               :: n_thermal_conductivity

    real(prec), dimension(nz) :: dry_cond
    real(prec), dimension(nz) :: sat_cond
! ----------------------------------------------------------------------
! Function: tfm_temperature_liquid_cond_Voigt
!
! The function computes the effective thermal conductivity of
! unsaturated firn following Kaviany (1991), using the model of
! Voigt (1889).
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of variables "density", "temperature", and
!     "liuqid_water".
!   density: Density along the firn profile (kg m**-3).
!   temperature: Temperature along the firn profile (K).
!   liquid_water: Liquid water content along the profile (1).
!   dry_thermcond_model: Model to use for the computation of effective
!     thermal conductvity at dry conditions.
!   sat_thermcond_model: Model to use for the computation of effective
!     thermal conductviity at saturated conditions.
! ----------------------------------------------------------------------

    dry_cond = dry_thermcond_model(nz, density, temperature)
    sat_cond = sat_thermcond_model(nz, density)

    n_thermal_conductivity = tfm_temperature_mixture_Voigt( &
    & nz, liquid_water, dry_cond, sat_cond                  &
    )
  end function tfm_temperature_liquid_cond_Voigt
end module tfm_temperature
