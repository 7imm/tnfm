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
  character(len=100) :: solve_liquid
  character(len=100) :: solve_grain_growth

  character(len=100) :: forcing_input_file
  character(len=100) :: initial_input_file

  real(prec) :: time_step
  real(prec) :: spinup

  real(prec), dimension(:), allocatable :: depth
  real(prec), dimension(:), allocatable :: density
  real(prec), dimension(:), allocatable :: temperature
  real(prec), dimension(:), allocatable :: grain_radius
  real(prec), dimension(:), allocatable :: age

  ! simulation properties
  real(prec), dimension(:,:), allocatable :: forcing
  real(prec), dimension(:), allocatable   :: runoff
  type(llProps)                           :: props

  real(prec), dimension(7) :: spinup_forcing

  type(sim_models) :: models

  integer    :: n
  integer    :: nt = 0
  integer    :: snt = 0
  integer    :: nz = 0
  integer    :: np = 0
  integer    :: tic, toc
  integer    :: t
  real(prec) :: rate


  ! configuration file list
  namelist /config/              &
  &  solve_density,              &
  &  solve_temperature,          &
  &  solve_heat_capacity,        &
  &  solve_thermal_conductivity, &
  &  solve_liquid,               &
  &  solve_grain_growth,         &
  &  forcing_input_file,         &
  &  initial_input_file,         &
  &  time_step,                  &
  &  spinup

  ! read configuration
  open(111, file=CONFIGURATION_FILE, action='read')
    read(unit=111, nml=config)
  close(111)

  ! sinup configuration
  if ( spinup > 0.0 ) then
    snt = int((spinup * SECONDS_YEAR) / time_step)
  else
    snt = 0
  end if

  ! forcing import
  call tfm_file_length(trim(forcing_input_file), nt)
  allocate(forcing(7,nt), runoff(nt))
  call tfm_read_csv(trim(forcing_input_file), 7, nt, forcing)

  ! init import
  call tfm_read_init(trim(initial_input_file), props)
  nz = props%depth%length
  np = nz + nt + snt

  ! model initialization
  call tfm_num_modelinit(                                         &
  &  solve_density=trim(solve_density),                           &
  &  solve_temperature=trim(solve_temperature),                   &
  &  solve_heat_capacity=trim(solve_heat_capacity),               &
  &  solve_thermal_conductivity=trim(solve_thermal_conductivity), &
  &  solve_liquid=trim(solve_liquid),                             &
  &  solve_grain_growth=trim(solve_grain_growth),                 &
  &  models=models                                                &
  )
  call llAppendData(props%heatcap, nz, models%heatcap_model(nz))
  call llAppendData(props%thermcond, nz, models%thermcond_model(nz, llGetData(props%density)))

  ! freedback
  print '(a,a)', 'Model Definitions'
  print '(a,a)', '================='
  print '(a,a)', 'density:                   ', trim(solve_density)
  print '(a,a)', 'temperature:               ', trim(solve_temperature)
  print '(a,a)', 'heat capacity:             ', trim(solve_heat_capacity)
  print '(a,a)', 'thermal conductivity:      ', trim(solve_thermal_conductivity)
  print '(a,a)', 'liquid:                    ', trim(solve_liquid)
  print '(a,a)', 'grain growth:              ', trim(solve_grain_growth)
  print '(a,a)', ''
  print '(a,a)', 'read forcing from:         ', trim(forcing_input_file)
  print '(a,a)', 'read initial profile from: ', trim(initial_input_file)
  print '(a,a)', ''
  print '(a,f13.2)', 'time step (s):            ', time_step
  print '(a,i10)',   'number of time steps:     ', nt
  print '(a,i10)',   'points in inital profile: ', nz
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
    &  (sum(forcing(6,:)) / nt)  & ! _> surface grain radius
    /)

    do t = 1, snt, 1

      call tfm_tools_indicate_tstep(snt, t)

      call tfm_num_surface( &
      &  np, time_step,     &
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

  ! time loop
  print '(a)', 'Simulation run:'
  do t = 1, nt, 1
    
    call tfm_tools_indicate_tstep(nt, t)

    call tfm_num_surface( &
    &  np, time_step,     &
    &  forcing(:,t),      &
    &  models,            &
    &  props             &
    )

    call tfm_num_step(         &
    &  time_step,              &
    &  models=models,          &
    &  props=props,            &
    &  runoff=runoff(t),       &
    &  liquid_acc=forcing(6,t) &
    )
  end do

  ! feedback
  call system_clock(toc, rate)
  print *, ''
  write(*, '(a,f10.2,a)') 'time elapsed: ', real(toc - tic) / real(rate), ' s'

  ! simple output
  nz = props%depth%length
  allocate(depth(nz), density(nz), temperature(nz), grain_radius(nz), age(nz))
  depth        = llGetData(props%depth)
  density      = llGetData(props%density)
  temperature  = llGetData(props%temperature)
  grain_radius = llGetData(props%grain_radius)
  age          = llGetData(props%age)
  open(333, file='tfm.out', status='replace', action='write')
  do n = nz, 1, -1
    write(333,*) depth(n), density(n), temperature(n), grain_radius(n), age(n)
  end do
  close(333)
  deallocate(depth, density, temperature, grain_radius, age)

  ! memoray deallocation
  call llFreeList(props%depth)
  call llFreeList(props%density)
  call llFreeList(props%temperature)
  call llFreeList(props%heatcap)
  call llFreeList(props%thermcond)
  call llFreeList(props%grain_radius)
  call llFreeList(props%liquidwater)
  call llFreeList(props%age)
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


subroutine tfm_read_init(input_file, props)
  use tfm_constants
  use tfm_llStructure
  implicit none

  character(len=*), intent(in) :: input_file
  type(llProps), intent(inout) :: props

  integer                  :: stat
  real(prec), dimension(6) :: init

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
end subroutine tfm_read_init
