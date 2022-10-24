module tfm_output
  use tfm_essentials
  use netcdf
  implicit none


  type out_ncid
    integer :: id
    integer :: time_dim
    integer :: depth_dim
    integer :: time_var
    integer :: depth_var
    integer :: dens_var
    integer :: temp_var
  end type out_ncid


  contains


  subroutine tfm_ncout_init(fname, deflate_level, np, ncid)
    implicit none

    character(len=*), intent(in)          :: fname
    integer, intent(in)                   :: deflate_level
    integer, intent(in)                   :: np
    type(out_ncid), intent(inout)         :: ncid
    integer                               :: stat

    stat = nf90_create(path=fname, cmode=nf90_netcdf4, ncid=ncid%id)

    stat = nf90_def_dim(ncid%id, 'time',  nf90_unlimited, ncid%time_dim)
    stat = nf90_def_dim(ncid%id, 'depth', np, ncid%depth_dim)

    stat = nf90_def_var(                    &
    &  ncid%id, 'time', nf90_float,         &
    &  (/ ncid%time_dim /),                 &
    &  ncid%time_var,                       &
    &  deflate_level=deflate_level          &
    )
    stat = nf90_def_var(                    &
    &  ncid%id, 'depth', nf90_float,        &
    &  (/ ncid%depth_dim, ncid%time_dim /), &
    &  ncid%depth_var,                      &
    &  deflate_level=deflate_level          &
    )
    stat = nf90_def_var(                    &
    &  ncid%id, 'density', nf90_float,      &
    &  (/ ncid%depth_dim, ncid%time_dim /), &
    &  ncid%dens_var,                       &
    &  deflate_level=deflate_level          &
    )
    stat = nf90_def_var(                    &
    &  ncid%id, 'temperature', nf90_float,  &
    &  (/ ncid%depth_dim, ncid%time_dim /), &
    &  ncid%temp_var,                       &
    &  deflate_level=deflate_level          &
    )

    stat = nf90_def_var_fill(ncid%id, ncid%time_var,  0, -9999.9)
    stat = nf90_def_var_fill(ncid%id, ncid%depth_var, 0, -9999.9)
    stat = nf90_def_var_fill(ncid%id, ncid%dens_var,  0, -9999.9)
    stat = nf90_def_var_fill(ncid%id, ncid%temp_var,  0, -9999.9)
  end subroutine tfm_ncout_init


  subroutine tfm_ncout_write(ncid, t, nz, props, time)
    implicit none

    type(out_ncid), intent(in)              :: ncid
    integer, intent(inout)                  :: t
    integer, intent(in)                     :: nz
    real(prec), dimension(8,nz), intent(in) :: props
    real(prec), intent(in)                  :: time
    
    integer :: stat

    t = t + 1

    stat = nf90_put_var(        &
    &  ncid%id, ncid%time_var,  &
    &  (/ time /),              &
    &  start=(/ t /),           &
    &  count=(/ 1 /)            &
    )
    stat = nf90_put_var(        &
    &  ncid%id, ncid%depth_var, &
    &  props(1,:),              &
    &  start=(/ 1, t /),        &
    &  count=(/ nz, 1 /)        &
    )
    stat = nf90_put_var(        &
    &  ncid%id, ncid%dens_var,  &
    &  props(2,:),              &
    &  start=(/ 1, t /),        &
    &  count=(/ nz, 1 /)        &
    )
    stat = nf90_put_var(        &
    &  ncid%id, ncid%temp_var,  &
    &  props(3,:),              &
    &  start=(/ 1, t /),        &
    &  count=(/ nz, 1 /)        &
    )
  end subroutine tfm_ncout_write


  subroutine tfm_ncout_close(ncid)
    implicit none

    type(out_ncid), intent(in) :: ncid
    integer                    :: stat

    stat = nf90_close(ncid%id)
  end subroutine tfm_ncout_close
end module tfm_output


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
