module tfm_temperature
  use tfm_essentials
  implicit none


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


  function tfm_temperature_capacity_paterson1994(nz) result(n_heat_capacity)
    implicit none
   
    integer, intent(in)       :: nz
    real(prec), dimension(nz) :: n_heat_capacity

    n_heat_capacity = 2009.0
  end function tfm_temperature_capacity_paterson1994


  function tfm_temperature_conduct_sturm2007(nz, density) &
    & result(n_thermal_conductivity)
    implicit none

    integer, intent(in)                   :: nz
    real(prec), dimension(nz), intent(in) :: density
    real(prec), dimension(nz)             :: n_thermal_conductivity

    n_thermal_conductivity = (         &
    &  (0.138)                         &
    &  - ((1.010e-3) * density)        &
    &  + ((3.233e-6) * (density**2.0)) &
    )
  end function tfm_temperature_conduct_sturm2007
end module tfm_temperature
