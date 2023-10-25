module tfm_liquid
  use tfm_constants
  use tfm_essentials 

  implicit none


  real(prec), parameter :: IMP_DENSITY = 830.0


  type vanGenuchtenParameters
    real(prec), dimension(:), allocatable :: alpha
    real(prec), dimension(:), allocatable :: n
    real(prec), dimension(:), allocatable :: m
  end type vanGenuchtenParameters

  
  ! interface for van Genuchten models
  interface
    subroutine van_genuchten_inter(nz, density, grain_radius, vg_params)
      use tfm_essentials
      import vanGenuchtenParameters
      implicit none

      integer, intent(in)                         :: nz
      real(prec), dimension(nz), intent(in)       :: density
      real(prec), dimension(nz), intent(in)       :: grain_radius
      type(vanGenuchtenParameters), intent(inout) :: vg_params
    end subroutine van_genuchten_inter
  end interface


  contains


  subroutine vgAllocateParams(nz, vg_params)
    implicit none

    integer, intent(in)          :: nz
    type(vanGenuchtenParameters) :: vg_params

    allocate(vg_params%alpha(nz))
    allocate(vg_params%n(nz))
    allocate(vg_params%m(nz))
  end subroutine vgAllocateParams


  subroutine vgDeallocateParams(vg_params)
    implicit none
    
    type(vanGenuchtenParameters) :: vg_params

    deallocate(vg_params%alpha)
    deallocate(vg_params%n)
    deallocate(vg_params%m)
  end subroutine vgDeallocateParams


  subroutine vgParametersYamaguchi2012(nz, density, grain_radius, vg_params)
    implicit none

    integer, intent(in)                         :: nz
    real(prec), dimension(nz), intent(in)       :: density
    real(prec), dimension(nz), intent(in)       :: grain_radius
    type(vanGenuchtenParameters), intent(inout) :: vg_params
! ----------------------------------------------------------------------
! Subroutine: vgParametersYamaguchi2012
!
! van Geuchten parameters for snow according to Yamaguchi et al. 2012.
!
! Yamaguchi, S., Watanabe, K., Katsushima, T., Sato, A., and Kumakura
! (2012). Dependence of the water retention curve of snow on snow
! characteristics. Annals of Glaciology, 53 (61), pp. 6-12,
! https://doi.org/10.3189/2012AoG61A001
!
! Author: Timm Schultz
! 
! Arguments:
!   nz: Dimension of variables "density", "grain radius", and
!     "vg_params".
!   density: Density along the profile (m/s).
!   grain_radius: Grain radius along the profile (kg/m**3).
!   vg_params - on input: van Genuchten parameters 
!     (of type vanGenuchtenParameters).
!
! Result:
!   vg_params - on output: van Genuchten parameters
!     (of type vanGenuchtenParameters).
! ----------------------------------------------------------------------

    ! alpha parameter
    vg_params%alpha = (                      &
    &  4.4E+6_prec                           &
    &  * ((                                  &
    &    density / (2.0_prec * grain_radius) &
    &  )**(-0.98_prec))                      &
    )

    ! n parameter
    vg_params%n = (                            &
    &  1.0_prec + (                            &
    &    2.7E-3_prec                           &
    &    * ((                                  &
    &      density / (2.0_prec * grain_radius) &
    &    )**(0.61_prec))                       &
    &  )                                       &
    )

    ! m parameter
    vg_params%m = (1.0_prec - (1.0_prec / vg_params%n))
  end subroutine vgParametersYamaguchi2012


  subroutine vgParametersDaanen2009(nz, density, grain_radius, vg_params)
    implicit none

    integer, intent(in)                         :: nz
    real(prec), dimension(nz), intent(in)       :: density
    real(prec), dimension(nz), intent(in)       :: grain_radius
    type(vanGenuchtenParameters), intent(inout) :: vg_params
! ----------------------------------------------------------------------
! Subroutine: vgParametersDaanen2009
!
! van Geuchten parameters for snow according to Daanen & Nieber 2009.
!
! Daanen, R. P. and Nieber, J. L. (2009). Model for Coupled Liquid Water
! Flow and Heat Transport with Phase Change in a Snowpack. Journal of
! Cold Regions Engineering, 23 (2), pp. 43-68,
! https://doi.org/10.1061/(ASCE)0887-381X(2009)23:2(43)
!
! Author: Timm Schultz
! 
! Arguments:
!   nz: Dimension of variables "density", "grain radius", and
!     "vg_params".
!   density: Density along the profile (m/s).
!   grain_radius: Grain radius along the profile (kg/m**3).
!   vg_params - on input: van Genuchten parameters 
!     (of type vanGenuchtenParameters).
!
! Result:
!   vg_params - on output: van Genuchten parameters
!     (of type vanGenuchtenParameters).
! ----------------------------------------------------------------------

    call tfm_essentials_do_nothing(nz, density)

    vg_params%alpha = (                                        &
    &  (30.0_prec * (2.0_prec * (grain_radius * 1000.0_prec))) &
    &  + 12.0_prec                                             &
    )

    vg_params%n = (                                            &
    &  (0.8_prec * (2.0_prec * (grain_radius * 1000.0_prec))) &
    &  + 3.0_prec                                              &
    )

    vg_params%m = (1.0_prec - (1.0_prec / vg_params%n))
  end subroutine vgParametersDaanen2009


  subroutine vgParametersYamaguchi2010(nz, density, grain_radius, vg_params)
    implicit none

    integer, intent(in)                         :: nz
    real(prec), dimension(nz), intent(in)       :: density
    real(prec), dimension(nz), intent(in)       :: grain_radius
    type(vanGenuchtenParameters), intent(inout) :: vg_params
! ----------------------------------------------------------------------
! Subroutine: vgParametersYamaguchi2010
!
! van Geuchten parameters for snow according to Yamaguchi et al. 2010.
!
! Yamaguchi, S., Katsushima, T., Sato, A. and Kumakura, T. (2010). Water
! retention curve of snow with different grain sizes. Cold Regions
! Science and Technology, 64 (2), pp. 87-93,
! https://doi.org/10.1016/j.coldregions.2010.05.008
!
! Author: Timm Schultz
! 
! Arguments:
!   nz: Dimension of variables "density", "grain radius", and
!     "vg_params".
!   density: Density along the profile (m/s).
!   grain_radius: Grain radius along the profile (kg/m**3).
!   vg_params - on input: van Genuchten parameters 
!     (of type vanGenuchtenParameters).
!
! Result:
!   vg_params - on output: van Genuchten parameters
!     (of type vanGenuchtenParameters).
! ----------------------------------------------------------------------

    call tfm_essentials_do_nothing(nz, density)

    vg_params%alpha = (                                  &
    &  (7.3_prec * (2.0 * (grain_radius * 1000.0_prec))) &
    &  + 1.9_prec                                        &
    )

    vg_params%n = (                                       &
    &  (-3.3_prec * (2.0 * (grain_radius * 1000.0_prec))) &
    &  + 14.4_prec                                        &
    )

    vg_params%m = (1.0_prec - (1.0_prec / vg_params%n))
  end subroutine vgParametersYamaguchi2010


  function vgDryLayers(nz, water_content) result(n_water_content)
    implicit none

    integer, intent(in)                   :: nz
    real(prec), dimension(nz), intent(in) :: water_content
    real(prec), dimension(nz)             :: n_water_content

    n_water_content = max(water_content, 1.0E-6_prec)
  end function vgDryLayers


  function vgSaturationWC(nz, density) result(saturation_wc)
    implicit none

    integer, intent(in)                   :: nz
    real(prec), dimension(nz), intent(in) :: density
    real(prec), dimension(nz)             :: saturation_wc

    saturation_wc = ((ICE_DENSITY - density) / WATER_DENSITY)
  end function vgSaturationWC


  function vgSaturationCondCalonne2012(nz, density, grain_radius) &
  
    & result(saturation_cond)
    implicit none

    integer, intent(in) :: nz
    real(prec), dimension(nz), intent(in) :: density
    real(prec), dimension(nz), intent(in) :: grain_radius
    real(prec), dimension(nz)             :: saturation_cond

    real(prec), parameter :: DYNAMIC_VISCOSITY = 0.001792_prec
! ----------------------------------------------------------------------
! Function: vgSaturationCondCalonne2012
!
! Hydraulic conductivity at saturation according to Calonne et al.
! (2012).
!
! Calonne, N., Geindreau, C., Flin, F., Morin, S., Lesaffre, B.,
! Rolland du Roscoat, S., and Charrier, P. (2012). 3-D image based
! numerical computations of snow permeability: links to specific surface
! area, density, and microstructual anisotropy. The Cryosphere, 6,
! pp. 939-951, https://doi.org/10.5194/tc-6-939-2012
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of variables "density" and "grain_radius".
!   density: Density along the firn profile (kg/m**3).
!   grain_raius: Equivalent grain radius along the profile (m).
!
! Result:
!  saturation_cond: Hydraulic conductivity at saturation (m/s).
! ----------------------------------------------------------------------

    saturation_cond = (                                    &
    &  ((WATER_DENSITY * ACC_GRAVItY) / DYNAMIC_VISCOSITY) &
    &  * (                                                 &
    &    0.75_prec                                         &
    &    * ((grain_radius / 1000.0_prec)**2.0_prec)        &
    &    * exp(-0.013_prec * (density))                    &
    &  )                                                   &
    )
  end function vgSaturationCondCalonne2012


  function vgSaturationCondShimizu1970(nz, density, grain_radius) &
    & result(saturation_cond)
    implicit none

    integer, intent(in)                   :: nz
    real(prec), dimension(nz), intent(in) :: density
    real(prec), dimension(nz), intent(in) :: grain_radius
    real(prec), dimension(nz)             :: saturation_cond

    real(prec), parameter :: DYNAMIC_VISCOSITY = 0.001792_prec
! ----------------------------------------------------------------------
! Function: vgSaturationCondShimizu1970
!
! Hydraulic conductivity at saturation according to Shimizu (1970).
!
! Shimizu, H. (1970). Air Permeability of Deposited Snow. Contributions
! from the Institute of Low Temperature Science, A22, pp. 1-32,
! http://hdl.handle.net/2115/20234
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of variables "density" and "grain_radius".
!   density: Density along the firn profile (kg/m**3).
!   grain_raius: Equivalent grain radius along the profile (m).
!
! Result:
!  saturation_cond: Hydraulic conductivity at saturation (m/s).
! ----------------------------------------------------------------------

    saturation_cond = (                                      &
    &  0.077_prec * (grain_radius**2.0_prec)                 &
    &  * exp(-0.0078_prec * density)                         &
    &  * ((WATER_DENSITY * ACC_GRAVITY) / DYNAMIC_VISCOSITY) &
    )
  end function vgSaturationCondShimizu1970


  function vgResidualWCWever2014(nz, water_content, residual_wc) &
    & result(n_residual_wc)
    implicit none

    integer, intent(in)                   :: nz
    real(prec), dimension(nz), intent(in) :: water_content
    real(prec), dimension(nz), intent(in) :: residual_wc
    real(prec), dimension(nz)             :: n_residual_wc
! ----------------------------------------------------------------------
! Function vgResidualWCWever2014
!
! Residual water content as defined by Wever et al. (2014).
!
! Wever, N., Fierz, N., Hirashima, H., and Lehnin, M. (2014). Solving
! Richards Equation for snow improves snowpack meltwater runoff
! estimations in detailed multi-layer snowpack model. The Cryosphere,
! 8, pp. 257-274, https://doi.org/10.5194/tc-8-257-2014
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of variables "water_content" and "residual_wc".
!   water_content: Volumetric water content along the profile (1).
!   residual_wc: Volumetric residual water content along the
!     profile (1).
!
! Result:
!   n_residual_wc: Volumetric residual water content along the
!     profile (1).
! ----------------------------------------------------------------------

    n_residual_wc = min(                             &
    &  0.02_prec,                                    &
    &  max((0.75_prec * water_content), residual_wc) &
    )
  end function vgResidualWCWever2014


  function vgRelativeHydraulicCond(nz, head, vg_params) &
    & result(rel_hydraulic_cond)
    implicit none

    integer, intent(in)                      :: nz
    real(prec), dimension(nz), intent(in)    :: head
    type(vanGenuchtenParameters), intent(in) :: vg_params
    real(prec), dimension(nz)                :: rel_hydraulic_cond

    rel_hydraulic_cond = (                                             &
    &  ((                                                              &
    &    1.0_prec                                                      &
    &    - (                                                           &
    &      ((vg_params%alpha * abs(head))**(vg_params%n - 1.0_prec))   &
    &      * ((                                                        &
    &        1.0_prec + ((vg_params%alpha * abs(head))**(vg_params%n)) &
    &      )**(-vg_params%m))                                          &
    &    )                                                             &
    &  )**2.0_prec)                                                    &
    &  / ((                                                            &
    &    1.0_prec + ((vg_params%alpha * abs(head))**(vg_params%n))     &
    &  )**(vg_params%m / 2.0_prec))                                    &
    )

    where ( head >= 0.0_prec )
      rel_hydraulic_cond = 1.0_prec
    end where
  end function vgRelativeHydraulicCond


  function vgSpecificMoistureCapNum(nz, head, residual_wc, saturation_wc, &
    & vg_params) result(specific_cap)
    implicit none

    integer, intent(in)                      :: nz
    real(prec), dimension(nz), intent(in)    :: head
    real(prec), dimension(nz), intent(in)    :: residual_wc
    real(prec), dimension(nz), intent(in)    :: saturation_wc
    type(vanGenuchtenParameters), intent(in) :: vg_params
    real(prec), dimension(nz)                :: specific_cap

    real(prec), parameter :: DHEAD = 1.0E-4_prec

    real(prec), dimension(nz) :: cplus
    real(prec), dimension(nz) :: cminus

    cplus = ( &
    &  (saturation_wc - residual_wc) &
    &  / (( &
    &    1.0_prec + ((vg_params%alpha * abs(head + DHEAD))**vg_params%n) &
    &  )**vg_params%m) &
    )
    cminus = ( &
    &  (saturation_wc - residual_wc) &
    &  / (( &
    &    1.0_prec + ((vg_params%alpha * abs(head - DHEAD))**vg_params%n) &
    &  )**vg_params%m) &
    )
    specific_cap = ((cplus - cminus) / (2.0 * DHEAD))

    where ( head >= 0.0_prec )
      specific_cap = 1.0_prec
    end where
  end function vgSpecificMoistureCapNum


  function vgWaterContent(nz, head, saturation_wc, residual_wc, vg_params) &
    & result(n_water_content)
    implicit none

    integer, intent(in)                      :: nz
    real(prec), dimension(nz), intent(in)    :: head
    real(prec), dimension(nz), intent(in)    :: saturation_wc
    real(prec), dimension(nz), intent(in)    :: residual_wc
    type(vanGenuchtenParameters), intent(in) :: vg_params
    real(prec), dimension(nz)                :: n_water_content

    n_water_content = (                                            &
    &  residual_wc                                                 &
    &  + (                                                         &
    &    (saturation_wc - residual_wc)                             &
    &    / ((                                                      &
    &      1.0_prec + ((vg_params%alpha * abs(head))**vg_params%n) &
    &    )**vg_params%m)                                           &
    &  )                                                           &
    )

    where ( head >= 0.0_prec )
      n_water_content = saturation_wc
    end where
  end function vgWaterContent


  function vgHead(nz, water_content, saturation_wc, residual_wc, &
    & vg_params) result(head)
    implicit none

    integer, intent(in)                      :: nz
    real(prec), dimension(nz), intent(in)    :: water_content
    real(prec), dimension(nz), intent(in)    :: saturation_wc
    real(prec), dimension(nz), intent(in)    :: residual_wc
    type(vanGenuchtenParameters), intent(in) :: vg_params
    real(prec), dimension(nz)                :: head

    head = (                                     &
    &  -(1.0_prec / vg_params%alpha)             &
    &  * ((                                      &
    &    ((                                      &
    &      (saturation_wc - residual_wc)         &
    &      / (water_content - residual_wc)       &
    &    )**(1.0_prec / vg_params%m)) - 1.0_prec &
    &  )**(1.0_prec / vg_params%n))              &
    )
  end function vgHead


  function vgLiquidMass(nz, depth, water_content) result(liquid_mass)
    implicit none

    integer, intent(in)                   :: nz
    real(prec), dimension(nz), intent(in) :: depth
    real(prec), dimension(nz), intent(in) :: water_content
    real(prec)                            :: liquid_mass

    liquid_mass = sum(               &
    &  (depth(2:nz) - depth(1:nz-1)) &
    &  * water_content(2:nz)         &
    &  * WATER_DENSITY               &
    )
  end function vgLiquidMass


  subroutine initKavetski2001(nz, &
    & water_content, saturation_wc, residual_wc, saturation_cond, vg_params, &
    & dt, truncerr_tolerance, &
    & head, dhdt, dwcdt, local_dt, safety_inp, eps_inp)
    implicit none
    
    integer, intent(in)                      :: nz
    real(prec), dimension(nz), intent(in)    :: water_content
    real(prec), dimension(nz), intent(in)    :: saturation_wc
    real(prec), dimension(nz), intent(in)    :: residual_wc
    real(prec), dimension(nz), intent(in)    :: saturation_cond
    type(vanGenuchtenParameters), intent(in) :: vg_params
    real(prec), intent(in)                   :: dt
    real(prec), dimension(2), intent(in)     :: truncerr_tolerance
    real(prec), intent(in), optional         :: safety_inp
    real(prec), intent(in), optional         :: eps_inp

    real(prec), dimension(nz), intent(inout) :: head
    real(prec), dimension(nz), intent(inout) :: dhdt
    real(prec), dimension(nz), intent(inout) :: dwcdt
    real(prec), intent(inout)                :: local_dt

    real(prec), dimension(nz) :: specific_cap
    real(prec), dimension(nz) :: hydraulic_cond
    real(prec)                :: safety
    real(prec)                :: eps

    if ( present(safety_inp) ) then
      safety = safety_inp
    else
      safety = 0.8_prec
    end if

    if ( present(eps_inp) ) then
      eps = eps_inp
    else
      eps = 1.0E-10_prec
    end if

    head = vgHead(nz, water_content, saturation_wc, residual_wc, vg_params)

    hydraulic_cond = (                                &
    &  saturation_cond                                &
    &  * vgRelativeHydraulicCond(nz, head, vg_params) &
    )

    specific_cap = vgSpecificMoistureCapNum(nz, head, residual_wc, &
      & saturation_wc, vg_params)

    dhdt = ((-hydraulic_cond * head) / specific_cap)
    dwcdt = (specific_cap * dhdt)

    local_dt = minval(safety * (                         &
    &  (                                                 &
    &    ((truncerr_tolerance(2)**0.5_prec) * abs(head)) &
    &    + (truncerr_tolerance(1)**0.5_prec)             &
    &  )                                                 &
    &  /(max(abs(dhdt), eps))                            &
    ))
    local_dt = min(dt, local_dt)
  end subroutine initKavetski2001


  subroutine timeStepKavetski2001(nz, local_dt, elapsed_time,      &
    & head, last_head, last_dhdt,                                  &
    & n_water_content, c_water_content, last_dwcdt,                &
    & truncerr_tolerance, safety_inp, rmin_inp, rmax_inp, eps_inp)
    implicit none

    integer, intent(in)                      :: nz
    real(prec), intent(inout)                :: local_dt
    real(prec), intent(inout)                :: elapsed_time
    real(prec), dimension(nz), intent(inout) :: head
    real(prec), dimension(nz), intent(inout) :: last_head
    real(prec), dimension(nz), intent(inout) :: last_dhdt
    real(prec), dimension(nz), intent(inout) :: n_water_content
    real(prec), dimension(nz), intent(inout) :: c_water_content
    real(prec), dimension(nz), intent(inout) :: last_dwcdt
    real(prec), dimension(2), intent(in)     :: truncerr_tolerance

    real(prec), intent(in), optional :: safety_inp
    real(prec), intent(in), optional :: rmin_inp
    real(prec), intent(in), optional :: rmax_inp
    real(prec), intent(in), optional :: eps_inp

    real(prec)                :: safety
    real(prec)                :: rmin
    real(prec)                :: rmax
    real(prec)                :: eps
    real(prec), dimension(nz) :: error
    real(prec), dimension(nz) :: condition
    real(prec)                :: dt_modifier
    integer                   :: icrit

    ! optional arguemnts
    if ( present(safety_inp) ) then
      safety = safety_inp
    else
      safety = 0.8_prec 
    end if

    if ( present(rmin_inp) ) then
      rmin = rmin_inp
    else
      rmin = 0.1_prec
    end if

    if ( present(rmax_inp) ) then
      rmax = rmax_inp
    else
      rmax = 2.0_prec
    end if

    if ( present(eps_inp) ) then
      eps = eps_inp
    else
      eps = 1.0E-10_prec
    end if

    ! time step adjustment
    error = (                                             &
    &  0.5_prec * local_dt                                &
    &  * abs(last_dhdt - ((head - last_head) / local_dt)) &
    )
    error(1) = 0.0_prec

    condition = (                            &
    &  error                                 &
    &  - (truncerr_tolerance(2) * abs(head)) &
    &  - (truncerr_tolerance(1))             &
    )

    icrit = maxloc(condition, 1)

    dt_modifier = safety * ((                                               &
    &  ((truncerr_tolerance(2) * abs(head(icrit))) + truncerr_tolerance(1)) &
    &  / (max(error(icrit), eps))                                           &
    )**0.5_prec)

    if ( condition(icrit) < 0.0_prec ) then

      last_dwcdt = (                                        &
      &  (                                                  &
      &    (2.0_prec * (c_water_content - n_water_content)) &
      &    - (local_dt * last_dwcdt)                        &
      &  )                                                  &
      &  / local_dt                                         &
      )
      last_dhdt = ((head - last_head) / local_dt)
      n_water_content = c_water_content
      c_water_content = n_water_content
      last_head = head
      elapsed_time = (elapsed_time + local_dt)
      local_dt = local_dt * min(dt_modifier, rmax)

    else if ( condition(icrit) >= 0.0_prec ) then

      head = last_head
      c_water_content = n_water_content
      local_dt = local_dt * max(dt_modifier, rmin)

    else
      print *, ''
      print *, '********************************************************'
      print *, '* Module: vgRichards                                   *'
      print *, '* Function: timeStepKavetski2011                       *'
      print *, '*                                                      *'
      print *, '* The mixed absolute-relative error, used to determine *'
      print *, '* whether the last time step is accepted or not, seems *'
      print *, '* to show an irregular value!                          *'
      print *, '* Stopping right here!                                 *'
      print *, '********************************************************'
      STOP
    end if
  end subroutine timeStepKavetski2001


  function solveRichardsEquation(nz, dt, depth, head, n_water_content,   &
    & c_water_content, last_dwcdt, hydraulic_cond, specific_cap, influx, &
    & outflux) result(d_head)
    implicit none

    integer, intent(in)                   :: nz
    real(prec), intent(in)                :: dt
    real(prec), dimension(nz), intent(in) :: depth
    real(prec), dimension(nz), intent(in) :: head
    real(prec), dimension(nz), intent(in) :: n_water_content
    real(prec), dimension(nz), intent(in) :: c_water_content
    real(prec), dimension(nz), intent(in) :: last_dwcdt
    real(prec), dimension(nz), intent(in) :: hydraulic_cond
    real(prec), dimension(nz), intent(in) :: specific_cap
    real(prec), intent(in)                :: influx
    real(prec), intent(in)                :: outflux
    real(prec), dimension(nz)             :: d_head

    real(prec), dimension(nz) :: alpha_u
    real(prec), dimension(nz) :: alpha_l
    real(prec), dimension(nz) :: am
    real(prec), dimension(nz) :: au
    real(prec), dimension(nz) :: al
    real(prec), dimension(nz) :: b

    integer    :: m
    real(prec) :: w
    real(prec) :: inter_cond

    alpha_u(2:nz-1) = (                                &
    &  (hydraulic_cond(2:nz-1) + hydraulic_cond(3:nz)) &
    &  / (                                             &
    &    (depth(3:nz) - depth(1:nz-2))                 &
    &    * (depth(3:nz) - depth(2:nz-1))               &
    &  )                                               &
    )
    alpha_u(1)  = -999999.9_prec
    alpha_u(nz) = -999999.9_prec

    alpha_l(2:nz-1) = (                                  &
    &  (hydraulic_cond(2:nz-1) + hydraulic_cond(1:nz-2)) &
    &  / (                                               &
    &    (depth(3:nz) - depth(1:nz-2))                   &
    &    * (depth(2:nz-1) - depth(1:nz-2))               &
    &  )                                                 &
    )
    alpha_l(1)  = -999999.9_prec
    alpha_l(nz) = -999999.9_prec

    ! main diagonal
    am(2:nz-1) = (                              &
    &  (2.0_prec * (specific_cap(2:nz-1)) / dt) &
    &  + alpha_u(2:nz-1)                        &
    &  + alpha_l(2:nz-1)                        &
    )

    ! upper diagonal
    au(2:nz-1) = (-alpha_u(2:nz-1))
    au(1)  = 0.0_prec
    au(nz) = -999999.9_prec

    ! lower diagonal
    al(2:nz-1) = (-alpha_l(2:nz-1))
    al(1)  = -999999.9_prec
    al(nz) = 0.0_prec

    ! right hand side
    b(2:nz-1) = (                                                              &
    &  ((2.0_prec * (n_water_content(2:nz-1) - c_water_content(2:nz-1))) / dt) &
    &  + (last_dwcdt(2:nz-1))                                                  &
    &  + (                                                                     &
    &    (hydraulic_cond(3:nz) - hydraulic_cond(1:nz-2))                       &
    &    / (depth(3:nz) - depth(1:nz-2))                                       &
    &  )                                                                       &
    &  - ((alpha_u(2:nz-1) + alpha_l(2:nz-1)) * head(2:nz-1))                  &
    &  - (-alpha_u(2:nz-1) * head(3:nz))                                       &
    &  - (-alpha_l(2:nz-1) * head(1:nz-2))                                     &
    )

    ! lower boundary condition (constant flux)
    if ( outflux >= 0.0_prec ) then

      inter_cond = (                                      &
      &  (hydraulic_cond(2) + hydraulic_cond(1))          &
      &  / (2.0_prec * ((depth(2) - depth(1))**2.0_prec)) &
      )

      au(1) = (-inter_cond)

      am(1) = (((2.0_prec * specific_cap(1)) / dt) + inter_cond)

      b(1) = (                                                           &
      &  - (+inter_cond * head(1))                                       &
      &  - (-inter_cond * head(2))                                       &
      &  + ((2.0_prec * (n_water_content(1) - c_water_content(1))) / dt) &
      &  + (last_dwcdt(1))                                               &
      &  + (inter_cond * (depth(2) - depth(1)))                          &
      &  + (outflux / (depth(2) - depth(1)))                             &
      )

    ! free surface
    else if ( outflux < 0.0_prec ) then

      am(1) = 1.0_prec
      au(1) = -1.0_prec
      b(1) = (head(2) - head(1))

    else
      print *, '*************************************************************'
      print *, '* Module: vgRichards                                        *'
      print *, '* Function: solveRichardsEquation                           *'
      print *, '*                                                           *'
      print *, '* The variables "outflux" seems to show an irregular value! *'
      print *, '* Stopping right here!                                      *'
      print *, '*************************************************************'
      STOP
    end if

    ! upper boundary condition (constant flux)
    if ( influx >= 0.0_prec ) then
      
      inter_cond = (                                          &
      &  (hydraulic_cond(nz) + hydraulic_cond(nz-1))          &
      &  / (2.0_prec * ((depth(nz) - depth(nz-1))**2.0_prec)) &
      )

      al(nz) = (-inter_cond)

      am(nz) = ((2.0_prec * (specific_cap(nz) / dt)) + inter_cond)
      
      b(nz) = (                                                            &
      &  - (-inter_cond * head(nz-1))                                      &
      &  - (+inter_cond * head(nz))                                        &
      &  + ((2.0_prec * (n_water_content(nz) - c_water_content(nz))) / dt) &
      &  + (last_dwcdt(nz))                                                &
      &  - (inter_cond * (depth(nz) - depth(nz-1)))                        &
      &  - (-influx / (depth(nz) - depth(nz-1)))                           &
      )

    ! free surface
    else if ( influx < 0.0_prec ) then

      am(nz) = 1.0_prec
      al(nz) = -1.0_prec
      b(nz) = (head(nz-1) - head(1))

    else
      print *, ''
      print *, '***********************************************************'
      print *, '* Module: vgRichards                                      *'
      print *, '* Function: solveRichardsEquation                         *'
      print *, '*                                                         *'
      print *, '* The variable "influx" seems to show an irregular value! *'
      print *, '* Stopping right here!                                    *'
      print *, '***********************************************************'
      STOP
    end if

    ! TDMA
    do m = 2, nz, 1
      w = al(m) / am(m-1)
      am(m) = am(m) - (w * au(m-1))
      b(m)  = b(m)  - (w * b(m-1))
    end do

    d_head(nz) = (b(nz) / am(nz))

    do m = (nz - 1), 1, -1
      d_head(m) = (b(m) - (au(m) * d_head(m+1))) / am(m)
    end do
  end function solveRichardsEquation


  function vgRichardsAdvanceTimeStep(nz, dt, depth, density, grain_radius, &
    & water_content, liquid_accumulation, van_genuchten_model)           &
    & result(n_water_content)
    implicit none

    integer, intent(in)                     :: nz
    real(prec), intent(in)                  :: dt
    real(prec), dimension(nz), intent(in)   :: depth
    real(prec), dimension(nz), intent(in)   :: density
    real(prec), dimension(nz), intent(in)   :: grain_radius
    real(prec), dimension(nz), intent(in)   :: water_content
    real(prec), intent(in)                  :: liquid_accumulation
    procedure(van_genuchten_inter), pointer :: van_genuchten_model
    real(prec), dimension(nz)               :: n_water_content

    integer    :: iter
    real(prec) :: backsteps
    real(prec) :: local_dt
    real(prec) :: elapsed_time
    real(prec) :: influx
    real(prec) :: outflux

    real(prec), dimension(2)  :: picard_tolerance
    real(prec), dimension(2)  :: step_tolerance
    real(prec), dimension(nz) :: residuum
    real(prec), dimension(nz) :: last_dhdt

    type(vanGenuchtenParameters) :: vg_params
    real(prec), dimension(nz)    :: saturation_wc
    real(prec), dimension(nz)    :: residual_wc
    real(prec), dimension(nz)    :: saturation_cond
    real(prec), dimension(nz)    :: rel_hydraulic_cond
    real(prec), dimension(nz)    :: hydraulic_cond
    real(prec), dimension(nz)    :: head
    real(prec), dimension(nz)    :: last_head
    real(prec), dimension(nz)    :: d_head
    real(prec), dimension(nz)    :: c_water_content
    real(prec), dimension(nz)    :: specific_cap
    real(prec), dimension(nz)    :: last_dwcdt
    real(prec), dimension(nz)    :: dry_layers
    real(prec), dimension(nz)    :: eff_saturation


    eff_saturation = 1.0E-3_prec
    saturation_wc = 0.9_prec * (1.0_prec - (density / ICE_DENSITY))

    where ( water_content == 0.0_prec )
      dry_layers = (                                             &
      &  (eff_saturation * saturation_wc)                        &
      &  / (1.0_prec - 0.75_prec + (0.75_prec * eff_saturation)) &
      )
    else where
      dry_layers = 0.0_prec
    end where

    ! defintions
    step_tolerance = (/ &
    &  0.0_prec,        &
    &  1.0E-2_prec      &
    /)
    picard_tolerance = (/             &
    &  0.01_prec * step_tolerance(1), &
    &  0.01_prec * step_tolerance(2)  &
    /)

    backsteps = 0.0_prec

    n_water_content = water_content + dry_layers
    influx = liquid_accumulation
    outflux = -1

    ! van Genuchten parameters
    call vgAllocateParams(nz, vg_params)
    call van_genuchten_model(nz, density, grain_radius, vg_params)

    ! props
    residual_wc = min(0.02_prec, (0.75_prec * n_water_content), saturation_wc)
    saturation_cond = vgSaturationCondShimizu1970(nz, density, grain_radius)

    call initKavetski2001(                                                      &
    &  nz,                                                                      &
    &  n_water_content, saturation_wc, residual_wc, saturation_cond, vg_params, &
    &  dt, step_tolerance, head, last_dhdt, last_dwcdt, local_dt                &
    )
    local_dt = 1.0E-6_prec
    last_head = head

    ! time loop
    elapsed_time = 0.0_prec
    do while ( elapsed_time /= dt )

      residuum = 999999.9_prec
      iter = 0

      ! Picard loop
      do while ( (maxval(residuum) >= 0.0_prec) .and. (iter < 100) )
        
        rel_hydraulic_cond = vgRelativeHydraulicCond(nz, head, vg_params)
        hydraulic_cond = (rel_hydraulic_cond * saturation_cond)

        specific_cap = vgSpecificMoistureCapNum(nz, head, residual_wc, &
          & saturation_wc, vg_params)

        d_head = solveRichardsEquation(nz, local_dt, depth, head,         &
          & n_water_content, c_water_content, last_dwcdt, hydraulic_cond, &
          & specific_cap, influx, outflux)

        head = (head + d_head)

        c_water_content = vgWaterContent(nz, head, &
          & saturation_wc, residual_wc, vg_params)

        ! exception head
        if ( ( any(isnan(head)) ) .or. ( any(abs(head) > huge(abs(head))) ) )then
          print *, ''
          print *, '**************************************************'
          print *, '* Module: tfm_liquid                             *'
          print *, '* Function: vgRichardsAdvanceTimeStep            *'
          print *, '*                                                *'
          print *, '* The variable "head" shows one or more NaN or   *'
          print *, '* Inifinity values! This should not be the case! *'
          print *, '* Stopping right here!                           *'
          print *, '**************************************************'
          STOP
        end if

        iter = (iter + 1)
        residuum = (                           &
        &  abs(d_head)                         &
        &  - (picard_tolerance(2) * abs(head)) &
        &  - (picard_tolerance(1))             &
        )
      end do

      ! time and time step control
      call timeStepKavetski2001(                       &
      &  nz, local_dt, elapsed_time,                   &
      &  head, last_head, last_dhdt,                   &
      &  n_water_content, c_water_content, last_dwcdt, &
      &  step_tolerance                                &
      )

      if ( (elapsed_time + local_dt) > dt ) then
        local_dt = (dt - elapsed_time)
      end if
    end do

    n_water_content = (n_water_content - dry_layers)
    where ( (n_water_content <= 1.0E-4) )
      n_water_content = 0.0_prec
    end where

    ! deallocation of the van Genuchten parameters
    call vgDeallocateParams(vg_params)
  end function vgRichardsAdvanceTimeStep


  subroutine tfm_liquid_RichardsEquation(nz, dt, depth, density, temperature, &
    & grain_radius, water_content, liquid_accumulation, runoff, van_genuchten_model)
    use tfm_constants
    implicit none

    integer, intent(in)                   :: nz
    real(prec), intent(in)                :: dt
    real(prec), dimension(nz), intent(in) :: depth
    real(prec), dimension(nz), intent(in) :: grain_radius
    real(prec), intent(in)                :: liquid_accumulation

    procedure(van_genuchten_inter), pointer :: van_genuchten_model

    real(prec), dimension(nz), intent(inout) :: density
    real(prec), dimension(nz), intent(inout) :: temperature
    real(prec), dimension(nz), intent(inout) :: water_content
    real(prec), intent(inout)                :: runoff

    real(prec), dimension(nz) :: dz
    real(prec), dimension(nz) :: water_mass
    real(prec), dimension(nz) :: refreeze_cap
    real(prec), dimension(nz) :: ice_cap
    real(prec), dimension(nz) :: storage
    real(prec), dimension(nz) :: d_temperature
    real(prec), dimension(nz) :: d_density


    if ( (liquid_accumulation > 0.0_prec) .or. any(water_content > 0.0_prec) ) then
      water_content = vgRichardsAdvanceTimeStep(nz, dt, depth, density, &
        & grain_radius, water_content, liquid_accumulation, van_genuchten_model)
    end if

    dz(1) = 1.0_prec
    dz(2:nz) = (depth(2:nz) - depth(1:nz-1))

    water_mass = (water_content * WATER_DENSITY * dz)

    if ( any(water_mass > 0.0_prec) ) then
      
      ! potential mass per square meter that might refreeze
      refreeze_cap = (                                                  &
      &  (SPECIFIC_HEAT_ICE * density * dz * (MELT_TEMP - temperature)) &
      &  / (LATENT_HEAT)                                                &
      )

      ! pore space in kg ice equivalent per square meter
      ice_cap = ((1.0_prec - (density / ICE_DENSITY)) * ICE_DENSITY * dz)

      storage = min(refreeze_cap, ice_cap, water_mass)

      ! temperature change
      d_temperature = (LATENT_HEAT / (SPECIFIC_HEAT_ICE * density * dz)) * storage
      !d_temperature(1) = 0.0_prec
      temperature = (temperature + d_temperature)

      ! density change
      d_density = (storage / dz)
      d_density(1) = 0.0_prec
      density = (density + d_density)

      water_content = water_content - (storage / WATER_DENSITY / dz)
      where ( water_content < 1.0E-10_prec )
        water_content = 0.0_prec
      end where
    end if
  end subroutine tfm_liquid_RichardsEquation


  subroutine tfm_liquid_bucket(nz, dt, depth, density, temperature, &
    & grain_radius, liquid_water, infiltration_rate, runoff, van_genuchten_model)
    implicit none

    integer, intent(in)                   :: nz
    real(prec), intent(in)                :: dt
    real(prec), dimension(nz), intent(in) :: depth
    real(prec), dimension(nz), intent(in) :: grain_radius
    real(prec), intent(in)                :: infiltration_rate

    procedure(van_genuchten_inter), pointer :: van_genuchten_model

    real(prec), dimension(nz), intent(inout) :: density
    real(prec), dimension(nz), intent(inout) :: temperature
    real(prec), dimension(nz), intent(inout) :: liquid_water
    real(prec), intent(inout)                :: runoff

    integer    :: n
    real(prec) :: water
    real(prec) :: dz
    real(prec) :: irr_water_content
    real(prec) :: ice_cap, refreeze_cap, imp_cap
    real(prec) :: storage

    call tfm_essentials_do_nothing(nz, grain_radius)

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
