program tfm_example
  use tfm_constants
  use tfm_essentials
  use tfm_num
  use tfm_tools
  use tfm_llStructure
  implicit none

  ! Parameters
  character(len=*), parameter :: CONFIGURATION_FILE = './tfm.conf'

  ! input variables read from configuration file
  character(len=100) :: solve_density
  character(len=100) :: solve_temperature
  character(len=100) :: solve_heat_capacity
  character(len=100) :: solve_thermal_conductivity
  character(len=100) :: solve_liquid_thermal_conductivity
  character(len=100) :: solve_saturation_thermal_conductivity
  character(len=100) :: solve_liquid
  character(len=100) :: solve_van_genuchten
  character(len=100) :: solve_grain_growth

  character(len=100) :: forcing_input_file
  character(len=100) :: initial_input_file

  real(prec) :: time_step
  real(prec) :: spinup
  real(prec) :: max_profile_length
  real(prec) :: max_profile_age

  type(sim_models)                        :: models
  type(llProps)                           :: props
  real(prec), dimension(:,:), allocatable :: forcing
  real(prec), dimension(:), allocatable   :: runoff
  real(prec), dimension(7)                :: spinup_forcing

  integer    :: nt = 0
  integer    :: snt = 0
  integer    :: tic, toc
  integer    :: t
  real(prec) :: rate


  ! configuration file list
  namelist /config/                            &
  &  solve_density,                            &
  &  solve_temperature,                        &
  &  solve_heat_capacity,                      &
  &  solve_thermal_conductivity,               &
  &  solve_liquid_thermal_conductivity,        &
  &  solve_saturation_thermal_conductivity,    &
  &  solve_liquid,                             &
  &  solve_van_genuchten,                      &
  &  solve_grain_growth,                       &
  &  forcing_input_file,                       &
  &  initial_input_file,                       &
  &  time_step,                                &
  &  spinup,                                   &
  &  max_profile_length,                       &
  &  max_profile_age

  ! read configuration
  open(111, file=CONFIGURATION_FILE, action='read')
    read(unit=111, nml=config)
  close(111)

  ! model initialization
  call tfm_num_modelinit(                                                               &
  &  solve_density=trim(solve_density),                                                 &
  &  solve_temperature=trim(solve_temperature),                                         &
  &  solve_heat_capacity=trim(solve_heat_capacity),                                     &
  &  solve_thermal_conductivity=trim(solve_thermal_conductivity),                       &
  &  solve_liquid_thermal_conductivity=trim(solve_liquid_thermal_conductivity),         &
  &  solve_saturation_thermal_conductivity=trim(solve_saturation_thermal_conductivity), &
  &  solve_liquid=trim(solve_liquid),                                                   &
  &  solve_van_genuchten=trim(solve_van_genuchten),                                     &
  &  solve_grain_growth=trim(solve_grain_growth),                                       &
  &  models=models                                                                      &
  )

  ! sinup configuration
  if ( spinup > 0.0 ) then
    snt = int((spinup * SECONDS_YEAR) / time_step)
  else
    snt = 0
  end if

  ! forcing and init import
  call tfm_file_length(trim(forcing_input_file), nt)
  allocate(forcing(7,nt), runoff(nt))
  call tfm_read_csv(trim(forcing_input_file), 7, nt, forcing)
  call tfm_read_init(trim(initial_input_file), props, models)

  ! freedback
  print '(a,a)', 'Model Definitions'
  print '(a,a)', '================='
  print '(a,a)', 'density:                         ', trim(solve_density)
  print '(a,a)', 'temperature:                     ', trim(solve_temperature)
  print '(a,a)', 'heat capacity:                   ', trim(solve_heat_capacity)
  print '(a,a)', 'thermal conductivity:            ', trim(solve_thermal_conductivity)
  print '(a,a)', 'liquid thermal conductivity:     ', trim(solve_liquid_thermal_conductivity)
  print '(a,a)', 'saturation thermal conductivity: ', trim(solve_saturation_thermal_conductivity)
  print '(a,a)', 'liquid:                          ', trim(solve_liquid)
  print '(a,a)', 'van Genuchten:                   ', trim(solve_van_genuchten)
  print '(a,a)', 'grain growth:                    ', trim(solve_grain_growth)
  print '(a,a)', ''
  print '(a,a)', 'read forcing from:         ', trim(forcing_input_file)
  print '(a,a)', 'read initial profile from: ', trim(initial_input_file)
  print '(a,a)', ''
  print '(a,f13.2)', 'time step (s):            ', time_step
  print '(a,i10)',   'number of time steps:     ', nt
  print '(a,i10)',   'points in inital profile: ', props%depth%length
  print *, ''

  ! optional spinup using mean values for forcing
  call system_clock(tic, rate)
  if ( spinup > 0.0 ) then

    print '(a)', 'Spin-Up:'
    
    ! mean forcing for spinupp
    spinup_forcing = (/          &
    &  0.0_prec,                 & ! -> world time
    &  0.0_prec,                 & ! -> model time
    &  (sum(forcing(3,:)) / nt), & ! -> surface temperature
    &  (sum(forcing(4,:)) / nt), & ! -> surface densiy
    &  (sum(forcing(5,:)) / nt), & ! -> surface accumulation (solid)
    &  0.0_prec,                 & ! -> surface accumulation (liquid)
    &  (sum(forcing(7,:)) / nt)  & ! -> surface grain radius
    /)

    do t = 1, snt, 1

      call tfm_tools_indicate_tstep(snt, t)

      if (max_profile_length > 0.0) then
        call tfm_num_trimProfileLength(props, max_profile_length)
      end if

      if (max_profile_age > 0.0) then
        call tfm_num_trimProfileAge(props, max_profile_age)
      end if


      call tfm_num_surface( &
      &  time_step,         &
      &  spinup_forcing,    &
      &  models,            &
      &  props              &
      )

      call tfm_num_step(              &
      &  time_step,                   &
      &  models=models,               &
      &  props=props,                 &
      &  runoff=runoff(1),            &
      &  liquid_acc=spinup_forcing(6) &
      )
    end do
    
    print '(a)', ''
  end if

  call simpleOutput(0, props)

  ! time loop
  print '(a)', 'Simulation run:'
  do t = 1, nt, 1
    
    call tfm_tools_indicate_tstep(nt, t)

    if (max_profile_length > 0.0) then
      call tfm_num_trimProfileLength(props, max_profile_length)
    end if

    if (max_profile_age > 0.0) then
      call tfm_num_trimProfileAge(props, max_profile_age)
    end if

    call tfm_num_surface( &
    &  time_step,         &
    &  forcing(:,t),      &
    &  models,            &
    &  props              &
    )

    call tfm_num_step(         &
    &  time_step,              &
    &  models=models,          &
    &  props=props,            &
    &  runoff=runoff(t),       &
    &  liquid_acc=forcing(6,t) &
    )
    
    call simpleOutput(t, props)
  end do

  ! feedback
  call system_clock(toc, rate)
  print *, ''
  write(*, '(a,f10.2,a)') 'time elapsed: ', real(toc - tic) / real(rate), ' s'

  ! memoray deallocation
  call llPropsFree(props)
  deallocate(forcing, runoff)
end program tfm_example


subroutine tfm_file_length(input_file, length)
  implicit none

  character(len=*), intent(in) :: input_file
  integer, intent(inout)       :: length

  integer :: stat

  length = 0

  open(111, file=input_file, action='read')
    do
      read(111, fmt='(a)', iostat=stat)
      if ( stat /= 0 ) EXIT
      length = length + 1
    end do
  close(111)
end subroutine tfm_file_length


subroutine tfm_read_csv(input_file, width, length, arr)
  use tfm_constants
  implicit none

  character(len=*), intent(in)                       :: input_file
  integer, intent(in)                                :: width
  integer, intent(in)                                :: length
  real(prec), dimension(width,length), intent(inout) :: arr 

  integer n

  open(111, file=input_file, action='read')
    do n = 1, length, 1
      read(111,*) arr(:,n)
    end do
  close(111)
end subroutine tfm_read_csv


subroutine tfm_read_init(input_file, props, models)
  use tfm_constants
  use tfm_num
  use tfm_llStructure
  implicit none

  character(len=*), intent(in) :: input_file
  type(llProps), intent(inout) :: props
  type(sim_models), intent(in) :: models

  integer                  :: stat
  integer                  :: nz
  real(prec), dimension(6) :: init

  real(prec), dimension(:), allocatable :: heatcap
  real(prec), dimension(:), allocatable :: thermcond

  open(111, file=input_file, action='read')
    do
      read(111, *, iostat=stat) init
      if ( stat /= 0 ) EXIT

      call llAppendData(props%depth,        1, init(1))
      call llAppendData(props%density,      1, init(2))
      call llAppendData(props%temperature,  1, init(3))
      call llAppendData(props%grain_radius, 1, init(4))
      call llAppendData(props%liquidwater,  1, init(5))
      call llAppendData(props%age,          1, init(6))
    end do
  close(111)

  nz = props%depth%length
  allocate(heatcap(nz), thermcond(nz))

  heatcap   = models%heatcap_model( &
  &  nz,                            &
  &  llGetData(props%density),      &
  &  llGetData(props%temperature),  &
  &  llGetData(props%liquidwater)   &
  )
  thermcond = models%thermcond_model( &
  &  nz,                              &
  &  llGetData(props%density),        &
  &  llGetData(props%temperature)     &
  )

  call llAppendData(props%heatcap, nz, heatcap)
  call llAppendData(props%thermcond, nz, thermcond)

  deallocate(heatcap, thermcond)
end subroutine tfm_read_init


subroutine simpleOutput(nt, props)
  use tfm_llStructure
  implicit none

  integer, intent(in)       :: nt
  type(llProps), intent(in) :: props
  character(len=6)          :: out_number
  character(len=24)         :: out_name

  integer :: n

  real(prec), dimension(props%depth%length)        :: depth
  real(prec), dimension(props%density%length)      :: density
  real(prec), dimension(props%temperature%length)  :: temperature
  real(prec), dimension(props%grain_radius%length) :: grain_radius
  real(prec), dimension(props%age%length)          :: age
  real(prec), dimension(props%liquidwater%length)  :: liquidwater

  depth        = llGetData(props%depth)
  density      = llGetData(props%density)
  temperature  = llGetData(props%temperature)
  grain_radius = llGetData(props%grain_radius)
  age          = llGetData(props%age)
  liquidwater  = llGetData(props%liquidwater)

  write (out_number, '(I0.6)') nt
  out_name = 'tfm_output/tfm'//out_number//'.out'

  open(333, file=out_name, status='replace', action='write')
    do n = props%depth%length, 1, -1
      write(333,*) depth(n), density(n), temperature(n), grain_radius(n), age(n), liquidwater(n)
    end do
  close(333)
end subroutine simpleOutput
