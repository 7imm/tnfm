module merra2
  use tfm_essentials
  use tfm_tools
  use netcdf
  
  implicit none
  integer, parameter :: NT = 24

  contains


  subroutine merra2_get_information(inp_dir, file_list, dims)
    implicit none

    character(len=*), intent(in)         :: inp_dir, file_list
    integer, dimension(3), intent(inout) :: dims
    integer                              :: stat, file_len = 0, ncid
    character(len=100)                   :: inp_file

    character(len=(len(inp_dir) + len(file_list) + 28)) :: command

    ! file list generation
    command = ( &
      & 'find ' //                   &
      & inp_dir //                   &
      & ' -name "*.nc" | sort > ' // &
      & file_list                    & 
    )
    call system(command)

    ! total dimension of time
    open(111, file='ncfiles.list', action='read')
    do
      read(111, fmt='(a)', iostat=stat)
      if ( stat /= 0 ) EXIT
      file_len = file_len + 1
    end do
    close(111)
    dims(1) = NT * file_len

    ! dimension of longitude and latitude
    open(111, file=file_list, action='read')
    read(111, '(a)') inp_file
    close(111)

    stat = nf90_open(inp_file, nf90_nowrite, ncid)
    call merra2_get_dim(ncid, 'lon', 3, dims(2))
    call merra2_get_dim(ncid, 'lat', 3, dims(3))
    stat = nf90_close(ncid)
  end subroutine merra2_get_information


  subroutine merra2_coordinates(file_list, longitude, latitude, lon_dim, lat_dim)
    implicit none

    character(len=*), intent(in)                  :: file_list
    integer, intent(in)                           :: lon_dim, lat_dim
    real(prec), dimension(lon_dim), intent(inout) :: longitude
    real(prec), dimension(lat_dim), intent(inout) :: latitude
    integer                                       :: stat, ncid
    character(len=100)                            :: inp_file

    ! first input file from list
    open(111, file=file_list, action='read')
    read(111, '(a)') inp_file
    close(111)

    ! read netcdf
    stat = nf90_open(inp_file, nf90_nowrite, ncid)
    call merra2_get_var1(ncid, 'lon', longitude, lon_dim)
    call merra2_get_var1(ncid, 'lat', latitude, lat_dim)
    stat = nf90_close(ncid)
  end subroutine merra2_coordinates


  subroutine merra2_coord_var(file_list, nkey, keys, lon, lat, dims, var)
    implicit none

    character(len=*), intent(in)                  :: file_list
    integer, intent(in)                           :: nkey
    character(len=*), dimension(nkey), intent(in) :: keys
    real(prec), intent(in)                        :: lon, lat
    integer, dimension(3), intent(in)             :: dims

    real(prec), dimension(nkey,dims(1)), intent(inout) :: var

    real(prec), dimension(dims(2)) :: longitude
    real(prec), dimension(dims(3)) :: latitude
    integer                        :: ix, iy
    character(len=100)             :: ncfile
    integer                        :: n, m, stat, ncid, k, l

    ! indices of the given coordinates
    call merra2_coordinates(file_list, longitude, latitude, &
      & dims(2), dims(3))
    call merra2_lon_lat_index(lon, lat, longitude, latitude, &
      & dims(2), dims(3), ix, iy)


    ! loop over and read all nc-files in dataset
    open(111, file=file_list, action='read')
    n = 1
    do
      ! exit if end of file list is reached
      read(111, fmt='(a)', iostat=stat) ncfile
      if ( stat /= 0 ) EXIT

      ! open nc-file
      stat = nf90_open(ncfile, nf90_nowrite, ncid)

      ! get all data from keys given by parameter "keys"
      k = (n - 1) * NT + 1
      l = (n + 0) * NT + 0
      do m = 1, nkey, 1
        ! special case of "julian" to get the julian day from time
        if ( keys(m) == 'julian' ) then
          call merra2_julian_day(ncid, 0.0_prec, var(m,k:l))
        ! all other cases handling MERRA2 variable names
        else
          call merra2_get_var3(ncid, keys(m), var(m,k:l), (/ ix, iy, NT /))
        end if
      end do

      stat = nf90_close(ncid)
      n = n + 1
    end do
    close(111)
  end subroutine merra2_coord_var


  !subroutine merra2_coord_input(file_list, lon, lat, dims, input)
  !  implicit none
  !  
  !  character(len=*), intent(in)                    :: file_list
  !  integer, dimension(3), intent(in)               :: dims
  !  real(prec), intent(in)                          :: lon, lat
  !  real(prec), dimension(9,dims(1)), intent(inout) :: input

  !  real(prec), dimension(dims(2)) :: longitude
  !  real(prec), dimension(dims(3)) :: latitude
  !  integer                        :: ix, iy, stat, ncid, k, l, n
  !  character(len=100)             :: ncfile

  !  ! indices of requested coordinates
  !  call merra2_coordinates(file_list, longitude, latitude, dims(2), dims(3))
  !  call merra2_lon_lat_index(lon, lat, longitude, latitude, dims(2), dims(3), ix, iy)

  !  ! read all nc-files in dataset
  !  open(111, file=file_list, action='read')
  !  n = 1
  !  do
  !    read(111, fmt='(a)', iostat=stat) ncfile
  !    if ( stat /= 0 ) EXIT

  !    stat = nf90_open(ncfile, nf90_nowrite, ncid)

  !    k = (n - 1) * NT + 1
  !    l = (n + 0) * NT + 0

  !    call merra2_julian_day(ncid, 0.0_prec, input(1,k:l))
  !    call merra2_get_var3(ncid, 'EVAP',     input(2,k:l), (/ ix, iy, NT /))
  !    call merra2_get_var3(ncid, 'PRECSNO',  input(3,k:l), (/ ix, iy, NT /))
  !    call merra2_get_var3(ncid, 'PRECTOT',  input(4,k:l), (/ ix, iy, NT /))
  !    call merra2_get_var3(ncid, 'QLML',     input(5,k:l), (/ ix, iy, NT /))
  !    call merra2_get_var3(ncid, 'SPEEDMAX', input(6,k:l), (/ ix, iy, NT /))
  !    call merra2_get_var3(ncid, 'TLML',     input(7,k:l), (/ ix, iy, NT /))
  !    call merra2_get_var3(ncid, 'TSH',      input(8,k:l), (/ ix, iy, NT /))
  !    call merra2_get_var3(ncid, 'VLML',     input(9,k:l), (/ ix, iy, NT /))

  !    stat = nf90_close(ncid)
  !    n = n + 1
  !  end do
  !  close(111)
  !  call system('rm ' // file_list)
  !end subroutine merra2_coord_input


  subroutine merra2_julian_day(ncid, offset, time)
    implicit none
    
    real(prec), intent(in) :: offset
    integer, intent(in)    :: ncid

    real(prec), dimension(NT), intent(inout) :: time

    real(prec)         :: stat, yyyy, mm, dd, jd
    character (len=10) :: date

    ! gregorian date
    stat = nf90_get_att(ncid, nf90_global, 'RangeBeginningDate', date)
    read(date(1:4), *) yyyy
    read(date(6:7), *) mm
    read(date(9:10), *) dd

    ! julian day
    call julian_day(yyyy, mm, dd, jd)
    call merra2_get_var1(ncid, 'time', time, NT)
    time = (((time / 60.0) - 11.5) / 24.0) + jd - offset
  end subroutine merra2_julian_day


  subroutine merra2_lon_lat_index(lon, lat, longitude, latitude, &
    & longitude_dim, latitude_dim, ix, iy)
    implicit none
      
    real(prec), intent(in)                           :: lon, lat
    real(prec), dimension(longitude_dim), intent(in) :: longitude
    real(prec), dimension(latitude_dim), intent(in)  :: latitude
    integer, intent(in)                              :: longitude_dim, latitude_dim

    integer, intent (inout) :: ix, iy

    ! find the longitude index
    do ix = 1, longitude_dim
      if (longitude(ix) == lon) EXIT
    end do

    ! find the latitude index
    do iy = 1, latitude_dim
       if (latitude(iy) == lat) EXIT
    end do
  end subroutine merra2_lon_lat_index


  subroutine merra2_get_dim(ncid, key, keylen, dimlen)
    implicit none
  
    integer, intent(in)           :: ncid, keylen
    character(keylen), intent(in) :: key

    integer, intent(inout) :: dimlen
  
    integer :: dimid, stat
  
    ! length of the dimension
    stat = nf90_inq_dimid(ncid, key, dimid)
    stat = nf90_inquire_dimension(ncid, dimid, len=dimlen)
  end subroutine merra2_get_dim
  

  subroutine merra2_get_var1(ncid, key, var, varlen)
     implicit none

     integer, intent(in)          :: ncid, varlen
     character(len=*), intent(in) :: key

     real(prec), dimension(varlen), intent(inout) :: var

     integer :: stat, varid

     stat = nf90_inq_varid(ncid, key, varid)
     stat = nf90_get_var(ncid, varid, var)
  end subroutine merra2_get_var1


  subroutine merra2_get_var3(ncid, key, var, varlen)
     implicit none

     integer, intent(in)               :: ncid
     integer, dimension(3), intent(in) :: varlen
     character(len=*), intent(in)      :: key

     real(prec), dimension(varlen(3)), intent(inout) :: var

     integer               :: stat, varid
     integer, dimension(3) :: s, c
     character(len=100)    :: units
     real(prec)            :: unit_factor
    
     s = (/ varlen(1), varlen(2), 1 /)
     c = (/ 1, 1, varlen(3) /)

     stat = nf90_inq_varid(ncid, key, varid)
     stat = nf90_get_var(ncid, varid, var, start=s, count=c)
     stat = nf90_get_att(ncid, varid, 'units', units)
     call merra2_unit_conversion(units, unit_factor)
     var = var * unit_factor
  end subroutine merra2_get_var3


  subroutine merra2_unit_conversion(units, factor)
    implicit none

    character(len=*), intent(in) :: units

    real(prec), intent(inout) :: factor

    factor = 1.0

    if ( units == 'kg m-2 s-1' ) factor = (1.0 / 1000.0)
    if ( units == '1' )          factor = 1.0
    if ( units == 'm s-1' )      factor = 1.0
    if ( units == 'K' )          factor = 1.0
  end subroutine merra2_unit_conversion
end module merra2
