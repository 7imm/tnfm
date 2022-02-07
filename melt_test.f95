program melt_test
  use settings
  use merra2
  use tfm_preprocessing
  use tools
  implicit none

  character(len=6), dimension(3), parameter :: keys = (/ &
  &  'julian',                                           &
  &  'TLML  ',                                           &
  &  'TSH   '                                            &
  /)
  character(len=12), parameter  :: file_list = 'ncfiles.list'
  character(len=100), parameter :: inp_dir = '/home/timm/Downloads/MERRA-2/'
  real(prec), parameter         :: dt = 3600.0

  integer                                     :: m ,n, x, y
  integer                                     :: len_ice
  integer, dimension(3)                       :: dims
  real(prec)                                  :: lon, lat
  real(prec), dimension(:), allocatable       :: longitude
  real(prec), dimension(:), allocatable       :: latitude
  real(prec), dimension(:,:), allocatable     :: input
  real(prec), dimension(:,:), allocatable   :: melt

  real(prec), dimension(:), allocatable :: lons, lats
  integer, dimension(:), allocatable    :: check, ind_ice


  call merra2_get_information(inp_dir, file_list, dims)

  allocate(                         &
  &  input(3,dims(1)),              &
  &  longitude(dims(2)),            &
  &  latitude(dims(3))              &
  )

  call merra2_coordinates(file_list, longitude, latitude, dims(2), dims(3))

  allocate(                     &
  &  lons(dims(2) * dims(3)),   &
  &  lats(dims(2) * dims(3)),   &
  &  check(dims(2) * dims(3)),  &
  )

  n = 1
  do x = 1, dims(2), 1
    do y = 1, dims(3), 1
      lons(n)  = longitude(x)
      lats(n)  = latitude(y)
      n = n + 1
    end do
  end do

  call ice_mask_check(                        &
  &  './GimpIceMask_90m_2015_v1.2_epsg4326.nc', &
  &  (dims(2) * dims(3)), lons, lats, check   &
  )

  len_ice = 0
  do n = 1, (dims(2) * dims(3)), 1
    if ( check(n) == 1 ) len_ice = len_ice + 1
  end do

  allocate(ind_ice(len_ice), melt(len_ice,dims(1)))

  m = 0
  do n = 1, (dims(2) * dims(3)), 1
    if ( check(n) == 1 ) then
      m = m + 1
      ind_ice(m) = n
    end if
  end do


  do n = 1, len_ice, 1
    print *, n
    lon = lons(ind_ice(n))
    lat = lats(ind_ice(n))
    call merra2_coord_var(file_list, 3, keys, lon, lat, dims, input)
    call tfm_pre_pdd(nt, dt, input(2,:), 'medley2020_greenland', melt(n,:))
  end do

  deallocate(input, melt, longitude, latitude, lons, lats, check, ind_ice)
end program melt_test
