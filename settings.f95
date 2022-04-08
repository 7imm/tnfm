module settings
  integer, parameter :: prec = selected_real_kind(8)


  type sim_props
    real(prec), dimension(:), pointer :: depth       => null()
    real(prec), dimension(:), pointer :: density     => null()
    real(prec), dimension(:), pointer :: temperature => null()
    real(prec), dimension(:), pointer :: heatcap     => null()
    real(prec), dimension(:), pointer :: thermcond   => null()
    real(prec), dimension(:), pointer :: liquidwater => null()
    real(prec), dimension(:), pointer :: age         => null()
  end type sim_props


  contains


  !subroutine warn_depth(error)
  !  implicit none
  !  integer, intent(inout) :: error

  !  print *, 'depth not given'
  !  error = 1
  !end subroutine warn_depth


  !subroutine warn_density(error)
  !  implicit none
  !  integer, intent(inout) :: error

  !  print *, 'density not given'
  !  error = 1
  !end subroutine warn_density


  !subroutine warn_temperature(error)
  !  implicit none
  !  integer, intent(inout) :: error

  !  print *, 'temperature not given'
  !  error = 1
  !end subroutine warn_temperature


  !subroutine warn_liquid_acc(error)
  !  implicit none
  !  integer, intent(inout) :: error

  !  print *, 'liquid accumulation not given'
  !  error = 1
  !end subroutine warn_liquid_acc


  !subroutine warn_runoff(error)
  !  implicit none
  !  integer, intent(inout) :: error

  !  print *, 'runoff not given'
  !  error = 1
  !end subroutine warn_runoff


  !subroutine warn_solid_acc(error)
  !  implicit none
  !  integer, intent(inout) :: error

  !  print *, 'solid accumulation not given'
  !  error = 1
  !end subroutine warn_solid_acc
end module settings
