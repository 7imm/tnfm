module tfm_tools
  use tfm_essentials 
  implicit none
  contains


  subroutine tfm_tools_indicate_tstep(nt, t)
    implicit none

    integer, intent(in) :: nt, t
    real                :: perc
    integer             :: n, n_perc

    ! some computations
    perc = (real(t) / real(nt)) * 100.0
    n_perc = floor(perc / 5.0)

    ! number of time steps
    write(*, '(a)', advance='no') '\rtimestep: '
    write(*, '(i6,a,i6)', advance='no') t, ' of ', nt

    ! bar
    write(*, '(a)', advance='no') ' |'
    do n = 1, 20, 1
       if ( n < n_perc ) then
         write(*, '(a)', advance='no') '='
       else if (n == n_perc) then
         write(*, '(a)', advance='no') '>'
       else
         write(*, '(a)', advance='no') ' '
       end if
    end do
    write(*, '(a)', advance='no') '|'

    ! percentenge
    write(*, '(f6.2,a)', advance='no') perc, ' %'

    ! done
    if ( t == nt ) then
      write(*, '(a)') ''
      write(*, '(a)') 'Done!'
    end if
  end subroutine tfm_tools_indicate_tstep


  subroutine julian_day(year, month, day, julian)
    implicit none
    
    real(prec), intent(in)    :: year, month, day
    real(prec), intent(inout) :: julian
    real(prec)                :: y, m, d, b

    if ( month > 2 ) then
      y = year
      m = month
    else
      y = year - 1.0
      m = month + 12.0
    end if

    d = day
    b = 2.0 - floor(y / 100.0) + floor(y / 400.0)

    julian = (                       &
      &   floor(365.25 * (y + 4716)) &
      & + floor(30.6001 * (m + 1.0)) &
      & + d + b - 1524.5             &
    )
  end subroutine julian_day


  subroutine gregorian_date(julian_day, year, month, day)
    implicit none
    
    real(prec), intent(in)    :: julian_day
    real(prec), intent(inout) :: year, month, day
    real(prec)                :: z, f, alpha, a, b, c, d, e

    z = floor(julian_day + 0.5)
    f = julian_day + 0.5 - z

    alpha = floor((z - 1867216.25) / 36524.25)

    a = z + 1.0 + alpha - (alpha / 4.0)
    b = a + 1524.0
    c = floor((b - 122.1) / 365.25)
    d = floor(365.25 * c)
    e = floor((b - d) / 30.6001)

    day = b - d - floor(30.6001 * e) + f

    if ( e <= 13 ) then
      month = e - 1.0
      year = c - 4716.0
    else
      month = e -13.0
      year = c - 4715.0
    end if
  end subroutine gregorian_date
end module tfm_tools
