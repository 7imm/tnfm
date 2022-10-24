module tfm_essentials
  integer, parameter :: prec = selected_real_kind(8)


  type sim_props
    real(prec), dimension(:), pointer :: depth        => null()
    real(prec), dimension(:), pointer :: density      => null()
    real(prec), dimension(:), pointer :: temperature  => null()
    real(prec), dimension(:), pointer :: heatcap      => null()
    real(prec), dimension(:), pointer :: thermcond    => null()
    real(prec), dimension(:), pointer :: liquidwater  => null()
    real(prec), dimension(:), pointer :: age          => null()
    real(prec), dimension(:), pointer :: grain_radius => null()
  end type sim_props


  contains


end module tfm_essentials
