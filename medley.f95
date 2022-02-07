module medley
  use settings
  use merra2
  use tfm_density
  use tfm_temperature
  use tfm_liquid
  use tfm_num
  use tfm_preprocessing
  use tfm_constants
  implicit none

  contains

  subroutine testing(inp_dir, lon, lat)
    implicit none

    integer, parameter               :: stride = 24 * 5
    character(len=12), parameter     :: file_list = 'ncfiles.list'
    character(len=8), dimension(9), parameter :: keys = (/ &
    &  'julian  ',                                         &
    &  'EVAP    ',                                         &
    &  'PRECSNO ',                                         &
    &  'PRECTOT ',                                         &
    &  'QLML    ',                                         &
    &  'SPEEDMAX',                                         &
    &  'TLML    ',                                         &
    &  'TSH     ',                                         &
    &  'VLML    '                                          &
    /)

    real(prec), intent(in)       :: lon, lat
    character(len=*), intent(in) :: inp_dir

    integer, dimension(3)                   :: dims, rci_info
    real(prec), dimension(:,:), allocatable :: input, forcing, downsampled
    real(prec), dimension(6)                :: rci_dates

    integer    :: nz = 201
    real(prec) :: init_z = -200.0

    integer    :: n, t, nt, np
    real(prec) :: dt
    real(prec) :: runoff

    integer    :: tic, toc
    real(prec) :: rate

    real(prec), dimension(:,:), allocatable :: props
    type(sim_models) :: models


    ! defintion of the RCI
    !rci_dates = (/ &
    !  & 1991.0, 1.0,   0.5 + (0.5 / 24.0), &
    !  & 1992.0, 12.0, 30.5 + (0.5 / 24.0)  &
    !/)
    rci_dates = (/ &
      & 1980.0, 1.0,   0.5 + (0.5 / 24.0), &
      & 1995.0, 12.0, 30.5 + (0.5 / 24.0)  &
    /)

    ! information about the dataset
    call merra2_get_information(inp_dir, file_list, dims)

    ! input import
    call system_clock(tic, rate)
    print *, "data import..."
    allocate(input(9,dims(1)))
    call merra2_coord_var(file_list, 9, keys, lon, lat, dims, input)
    call system_clock(toc)
    print *, real(toc - tic) / real(rate)

    ! forcing generation 
    call system_clock(tic, rate)
    print *, "forcing generation..."
    allocate(forcing(6,dims(1)))
    call medley_forcing(input, dims(1), forcing)
    deallocate(input)

    ! downsampling
    print *, "downsampling..."
    nt = (dims(1) / stride)
    allocate(downsampled(6,nt))
    call medley_downsampling(stride, dims(1), forcing, nt, downsampled)
    deallocate(forcing)
    dt = downsampled(2,2) - downsampled(2,1)

    ! spin up information
    print *, "spin up information..."
    call medley_spinup(downsampled, nt, rci_dates, rci_info)

    ! model initialization
    print *, "model initialization..."
    np = ((rci_info(2) - rci_info(1)) * rci_info(3)) + nt + nz
    allocate(props(6,np))

    call tfm_num_modelinit(                    &
    &  solve_density='medley2020',             &
    &  solve_temperature='true',               &
    &  solve_heat_capacity='paterson1994',     &
    &  solve_thermal_conductivity='sturm2007', &
    &  solve_liquid='bucket',                  &
    &  models=models                           &
    )
    call medley_model_initialization(nt, np, nz, init_z, &
      & downsampled, props, models)

    ! spin up
    print *, "spin up..."
    do n = 1, rci_info(3), 1
      do t = rci_info(1), rci_info(2)
        call tfm_num_surface(nz, np, dt, downsampled(:,t), models, props)

        call tfm_num_step(              &
        &  nz, dt,                      &
        &  models=models,               &
        &  props=props(:,1:nz),         &
        &  runoff=runoff,               &
        &  liquid_acc=downsampled(6,t), &
        &  solid_acc=downsampled(5,t)   &
        )
      end do
    end do

    ! transient run
    props(1,1:nz) = props(1,1:nz) - props(1,nz)
    print *, "transient run..."
    do t = 1, nt, 1
      call tfm_num_surface(nz, np, dt, downsampled(:,t), models, props)

      call tfm_num_step(              &
      &  nz, dt,                      &
      &  models=models,               &
      &  props=props(:,1:nz),         &
      &  runoff=runoff,               &
      &  liquid_acc=downsampled(6,t), &
      &  solid_acc=downsampled(5,t)   &
      )
    end do
    call system_clock(toc)
    print *, real(toc - tic) / real(rate)

    open(333, file='depth_test.out', action='write')
    write (333,*) props(1,:)
    close(333)

    open(444, file='density_test.out', action='write')
    write (444,*) props(2,:)
    close(444)

    open(555, file='temperature_test.out', action='write')
    write (555,*) props(3,:)
    close(555)

    ! memory deallocation
    deallocate(downsampled, props)
  end subroutine testing 


  subroutine medley_model_initialization(nt, np, nz, init_z, &
    & forcing, props, models)
    implicit none

    integer, intent(in)                     :: nt, nz, np
    real(prec), intent(in)                  :: init_z
    real(prec), dimension(6,nt), intent(in) :: forcing
    type(sim_models), intent(in)            :: models

    real(prec), dimension(6,np), intent(inout) :: props

    integer    :: n
    real(prec) :: mean_temp, mean_surf_dens, mean_acc

    ! init
    props = -9999.9

    mean_temp      = sum(forcing(3,:)) / nt
    mean_surf_dens = sum(forcing(4,:)) / nt
    mean_acc       = sum(forcing(5,:) + forcing(6,:)) / nt

    ! initial profile (Herron & Langway)
    do n = 1, nz, 1
      props(1,n) = (nz - n) * (init_z / (nz - 1))
    end do

    call tfm_density_herron1980_analytical(nz, mean_temp, mean_acc, &
      & mean_surf_dens, props(1,1:nz), props(2,1:nz))

    props(3,1:nz) = mean_temp
    props(4,1:nz) = models%heatcap_model(nz)
    props(5,1:nz) = models%thermcond_model(nz, props(2,1:nz))
    props(6,1:nz) = 0.0
  end subroutine medley_model_initialization


  subroutine medley_downsampling(stride, nt_old, forcing_old, nt_new, forcing_new)
    implicit none

    integer, intent(in)                            :: stride, nt_old, nt_new
    real(prec), dimension(6,nt_old), intent(in)    :: forcing_old
    real(prec), dimension(6,nt_new), intent(inout) :: forcing_new

    integer    :: n, m, k, l, d

    forcing_new(1,1) = forcing_old(1,1)
    forcing_new(2,1) = forcing_old(2,1)
    forcing_new(3,1) = sum(forcing_old(3,1:(stride / 2))) / (stride / 2)
    forcing_new(4,1) = sum(forcing_old(4,1:(stride / 2))) / (stride / 2)
    forcing_new(5,1) = (sum(forcing_old(5,1:(stride / 2))) * 2.0) / stride
    forcing_new(6,1) = (sum(forcing_old(6,1:(stride / 2))) * 2.0) / stride

    do n = 2, nt_new
      d = ((n - 1) * stride) + 1
      k = d - (stride / 2)
      l = d + (stride / 2)

      do m = 1, 2
        forcing_new(m,n) = forcing_old(m,d)
      end do
      do m = 3, 6
        forcing_new(m,n) = sum(forcing_old(m,k:l)) / stride
      end do
    end do
  end subroutine medley_downsampling


  subroutine medley_forcing(input, nt, forcing)
    implicit none

    integer, intent(in)                        :: nt
    real(prec), dimension(9,nt), intent(in)    :: input
    real(prec), dimension(6,nt), intent(inout) :: forcing

    real(prec)                :: dt
    real(prec), dimension(nt) :: melt
    real(prec), dimension(nt) :: rain
    real(prec), dimension(nt) :: snowfall
    real(prec), dimension(nt) :: evaporation
    real(prec), dimension(nt) :: surfdensity

    ! temporary variables
    dt = (input(2,2) - input(2,1))

    call tfm_pre_pdd(nt, dt, input(7,:), 'medley2020_greenland', melt)
    call medley_forcing_surfdensity(input, nt, surfdensity)
 
    rain        = input(4,:) - input(3,:)
    snowfall    = input(3,:)
    evaporation = input(2,:)

    ! forcing array
    forcing(1,:) = input(1,:)                              ! date
    forcing(2,:) = (input(1,:) - input(1,1)) * SECONDS_DAY ! model time (seconds)
    forcing(3,:) = input(8,:)                              ! temperature
    forcing(4,:) = surfdensity                             ! surface density
    forcing(5,:) = snowfall - evaporation - melt           ! solid accumulation
    forcing(6,:) = melt + rain                             ! liquid accumulation
  end subroutine medley_forcing


  subroutine medley_forcing_surfdensity(input, nt, surf_dens)
    implicit none
      
    integer, intent(in)                     :: nt
    real(prec), dimension(9,nt), intent(in) :: input

    real(prec), dimension(nt), intent(inout) :: surf_dens

    real(prec) :: wind_north_mean, wind_max_mean, humidity_mean
    real(prec) :: accumulation_mean, temperature_mean

    ! mean values of the input data
    wind_north_mean  = (sum(input(9,:)) / nt)
    wind_max_mean    = (sum(input(6,:)) / nt)
    humidity_mean    = (sum(input(5,:)) / nt)
    temperature_mean = (sum(input(8,:)) / nt)

    ! accumulation is total precipitation - evaporation (- runoff)
    accumulation_mean = (sum(input(4,:) - input(2,:)) / nt) / WATER_DENSITY

    call tfm_pre_surfdens_medley2020(wind_north_mean, wind_max_mean, &
      & humidity_mean, accumulation_mean, temperature_mean, surf_dens(1))
    surf_dens(:) = surf_dens(1)
  end subroutine medley_forcing_surfdensity 


  subroutine medley_spinup(forcing, nt, rci_dates, rci_info)
    implicit none
    
    integer, intent(in)                     :: nt
    real(prec), dimension(6,nt), intent(in) :: forcing
    real(prec), dimension(6), intent(in)    :: rci_dates

    integer, dimension(3), intent(inout) :: rci_info

    real(prec) :: rci_start_jdate, rci_end_jdate
    integer    :: rci_start, rci_end

    ! defintion of the RCI interval
    call julian_day(rci_dates(1), rci_dates(2), rci_dates(3), rci_start_jdate)
    call julian_day(rci_dates(4), rci_dates(5), rci_dates(6), rci_end_jdate)

    do rci_start = 1, nt
      if ( forcing(1,rci_start) >= rci_start_jdate ) EXIT
    end do
    rci_info(1) = rci_start

    do rci_end = 1, nt
      if ( forcing(1,rci_end) >= rci_end_jdate) EXIT
    end do
    rci_info(2) = rci_end - 1

    ! Herron & Langway (1980) estimation
    call medley_herronlangway(forcing, nt, rci_start, rci_end, rci_info(3))
  end subroutine medley_spinup


  subroutine medley_herronlangway(forcing, nt, rci_start, rci_end, rci_cycles)
    implicit none

    integer, parameter    :: NZ = 1000
    real(prec), parameter :: MAXZ = -200.0
    
    integer, intent(in)                     :: nt, rci_start, rci_end
    real(prec), dimension(6,nt), intent(in) :: forcing

    integer, intent(inout) :: rci_cycles

    integer                   :: n
    real(prec)                :: surf_dens_mean, temp_mean, acc_mean
    real(prec)                :: burial_rate
    real(prec), dimension(nz) :: depth, density
    
    ! mean values over the RCI
    temp_mean      = (sum(forcing(3,rci_start:rci_end)) / (rci_end - rci_start))
    surf_dens_mean = (sum(forcing(4,rci_start:rci_end)) / (rci_end - rci_start))
    burial_rate    = (sum(forcing(5,rci_start:rci_end)) / (rci_end - rci_start))
    acc_mean       = (                                                  &
    &  sum(forcing(5,rci_start:rci_end) + forcing(6,rci_start:rci_end)) &
    &  / (rci_end - rci_start)                                          &
    )

    ! depth
    do n = 1, NZ, 1
      depth(n) = (nz - n) * (maxz / (nz - 1))
    end do

    call tfm_density_herron1980_analytical(NZ, temp_mean, acc_mean, surf_dens_mean, &
      & depth, density)

    ! depth of 910 kg m**-3
    do n = NZ, 1, -1
      if ( density(n) > 910.0 ) EXIT
    end do

    ! cycles of the rci needed to match the spin up
    rci_cycles = ceiling(                            &
    &    ((-depth(n) / burial_rate))                 &
    &  / (forcing(2,rci_end) - forcing(2,rci_start)) &
    )
  end subroutine medley_herronlangway
end module medley
