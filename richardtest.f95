program richardtest
  use settings
  use tfm_liquid
  implicit none

  real(prec) :: density
  real(prec) :: grain_radius
  real(prec) :: water_content

  density = 400.0
  grain_radius = 0.005
  water_content = 0.35

  call tfm_liquid_conduct(density, grain_radius, water_content)
end program richardtest
