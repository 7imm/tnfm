module tfm_tools
  use settings
  implicit none
  contains


  subroutine ice_mask_check(ice_mask_nc, nxy, lons, lats, check)
    use netcdf
    implicit none
    
    character(len=*), intent(in)           :: ice_mask_nc
    integer, intent(in)                    :: nxy
    real(prec), dimension(nxy), intent(in) :: lons, lats
    integer, dimension(nxy), intent(inout) :: check

    integer    :: stat, n
    integer    :: ncid, xid, yid, bid
    integer    :: xdim, ydim
    integer    :: ix, iy
    real(prec) :: neighbor

    real(prec), dimension(:), allocatable  :: longitude
    real(prec), dimension(:), allocatable  :: latitude
    real(prec), dimension(:), allocatable  :: xdiff, ydiff
    integer*8, dimension(:,:), allocatable :: band1

    stat = nf90_open(ice_mask_nc, nf90_nowrite, ncid)

    ! dimensions
    stat = nf90_inq_dimid(ncid, 'lon', xid)
    stat = nf90_inquire_dimension(ncid, xid, len=xdim)

    stat = nf90_inq_dimid(ncid, 'lat', yid)
    stat = nf90_inquire_dimension(ncid, yid, len=ydim)

    ! memory allocation
    allocate(longitude(xdim), latitude(ydim), band1(ydim,xdim), &
      & xdiff(xdim), ydiff(ydim))

    ! coordinates
    stat = nf90_inq_varid(ncid, 'lon', xid)
    stat = nf90_get_var(ncid, xid, longitude)

    stat = nf90_inq_varid(ncid, 'lat', yid)
    stat = nf90_get_var(ncid, yid, latitude)

    stat = nf90_inq_varid(ncid, 'Band1', bid)
    stat = nf90_get_var(ncid, bid, band1, start=(/ 1, 1 /), count=(/ xdim, ydim /))

    ! find nearest coordinate
    check = 0
    do n = 1, nxy, 1

      xdiff = abs(longitude - lons(n))
      neighbor = minval(xdiff)
      do ix = 1, xdim, 1
        if ( xdiff(ix) == neighbor ) EXIT
      end do

      ydiff = abs(latitude - lats(n))
      neighbor = minval(ydiff)
      do iy = 1, ydim, 1
        if ( ydiff(iy) == neighbor ) EXIT
      end do

      !check(n) = band1(iy,ix)
      if ( band1(iy,ix) == 1 ) check(n) = 1
    end do

    ! memory deallocation
    stat = nf90_close(ncid)
    deallocate(longitude, latitude, band1, xdiff, ydiff)
  end subroutine ice_mask_check


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
