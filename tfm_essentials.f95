module tfm_essentials
  use tfm_constants
  implicit none

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
! Subroutne: tfm_essentials_do_nothing
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
  

  subroutine tfm_essentials_mean_acc(nz, depth, density, age, mean_acc)
    use tfm_constants
    implicit none

    integer, intent(in)                      :: nz
    real(prec), dimension(nz), intent(in)    :: depth
    real(prec), dimension(nz), intent(in)    :: density
    real(prec), dimension(nz), intent(in)    :: age
    real(prec), dimension(nz), intent(inout) :: mean_acc

    integer :: n

!-----------------------------------------------------------------------
! Subroutine : tfm_essentials_mean_acc
!
! The subroutine computes the mean accumulation rate from the age of
! the firn profile. The concept follow the idea of calculating the
! mean accumulation rate over the life time of a firn parcel.
!
! See for example:
! Stevens, C. M., Verjans, V., Luding, J. M. D., Kahle, E. C.,
! Horlings, A. N., Horlings, B. I., and Waddington, E. D. The Community
! Firn Model (CFM) v1.0. Geosci. Model. Dev., 13, 4355-4377, (2020).
! https://doi.org/10.5194/gmd-13-4355-2020
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of variables "depth", "density", "age", "mean_acc".
!   depth: Depth of the firn profile (m).
!   density: Density of the firn profile (kg m**-3).
!   age: Age of the firn profile (s).
!   mean_acc - on input: Variable to store the mean accumulation rate.
!
! Result:
!   mean_acc - on output: Mean accumulation rate (m weq. a**-1).
!-----------------------------------------------------------------------
    
    mean_acc(nz) = 0.0
    do n = nz - 1, 1, -1
      mean_acc(n) = (depth(n+1) - depth(n)) * (density(n) / WATER_DENSITY)
      mean_acc(n) = mean_acc(n) + mean_acc(n+1)
    end do
  
    mean_acc(1:nz-1) = mean_acc(1:nz-1) / (age(1:nz-1) / SECONDS_YEAR)
  end subroutine tfm_essentials_mean_acc
end module tfm_essentials
