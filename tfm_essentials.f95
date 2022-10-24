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

  subroutine tfm_essentials_do_nothing(nz, variable)
    implicit none

    integer, intent(in)                   :: nz
    real(prec), dimension(nz), intent(in) :: variable
    real(prec), dimension(nz)             :: nothing
!-----------------------------------------------------------------------
! subroutne: tfm_essentials_do_nothing
!
! Routine to do "nothing" with a given variable. Used to avoid warning
! at compile time.
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of the given variable.
!   variable: Variable to do nothing about.
!-----------------------------------------------------------------------
    
    nothing = variable
  end subroutine tfm_essentials_do_nothing
end module tfm_essentials
