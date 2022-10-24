program kohnen03
  use tfm_essentials
  use tfm_num
  use tfm_tools
  implicit none

  character(len=*), parameter :: INIT_INP = './kohnen04init.dat'
  character(len=*), parameter :: FORC_INP = './kohnen04longForcing.dat'
  real(prec), parameter       :: DT = 3600.0_prec * 24.0_prec * (365.0_prec / 12.0_prec)

  real(prec), dimension(:,:), allocatable :: forcing
  real(prec), dimension(:,:), allocatable :: props
  real(prec), dimension(:), allocatable   :: runoff

  type(sim_models) :: models

  real(prec) :: rate
  
  integer :: tic, toc
  integer :: n
  integer :: np = 0
  integer :: nz = 0
  integer :: nt = 0
  integer :: t = 0


  ! model initialization
  call tfm_num_modelinit(                    &
  &  solve_density='arthern2010',            &
  &  solve_temperature='true',               &
  &  solve_heat_capacity='paterson1994',     &
  &  solve_thermal_conductivity='sturm2007', &
  &  solve_liquid='bucket',                  &
  &  solve_grain_growth='li2002',            &
  &  models=models                           &
  )

  ! gather information
  call kohnen03_sim_info(INIT_INP, FORC_INP, nz, nt)
  np = nz + nt

  ! forcing import
  allocate(forcing(7,nt), runoff(nt))
  call kohnen03_import_forcing(FORC_INP, nt, forcing)
  forcing(1,:) = 0.0_prec
  forcing(6,:) = 0.0_prec
  forcing(7,:) = 0.0005_prec
  
  ! init import
  allocate(props(8,np))
  props = -999999.9
  call kohnen03_import_init(INIT_INP, nz, props(:,1:nz))
  props(6,:) = 0.0_prec
  props(8,:) = 0.0005_prec

  ! output file
  open(333, file='kohnen04.out', status='replace', action='write')
  !do n = 1, nz, 1
  !  write(333,*) props(1,n), props(2,n), props(3,n)
  !end do
  !write(333,*) '###'

  ! time loop
  call system_clock(tic, rate)
  do t = 1, nt, 1
    
    call tfm_tools_indicate_tstep(nt, t)

    call tfm_num_surface( &
    &  nz, np, DT,        &
    &  forcing(:,t),      &
    &  models,            &
    &  props              &
    )

    call tfm_num_step(         &
    &  nz, dt,                 &
    &  models=models,          &
    &  props=props(:,1:nz),    &
    &  runoff=runoff(t),       &
    &  liquid_acc=forcing(6,t) &
    )

    !if ( t == 7000 ) EXIT

  end do

  ! simple output
  do n = 1, nz, 1
    write(333,*) props(1,n), props(2,n), props(3,n)
  end do
  write(333,*) '###'

  ! output file
  close(333)

  ! feedback
  call system_clock(toc, rate)
  print *, ''
  write(*, '(a,f10.2,a)') 'time elapsed: ', real(toc - tic) / real(rate), ' s'

  ! memory deallocation
  deallocate(forcing, props, runoff)
end program kohnen03


subroutine kohnen03_import_init(init_file, nz, props)
  use tfm_essentials
  implicit none

  character(len=*), intent(in)               :: init_file
  integer, intent(in)                        :: nz
  real(prec), dimension(8,nz), intent(inout) :: props

  integer n

  open(111, file=init_file, action='read')
    do n = 1, nz, 1
      read(111,*) props(1,n), props(2,n), props(3,n), props(7,n)
    end do
  close(111)
end subroutine kohnen03_import_init


subroutine kohnen03_import_forcing(forc_file, nt, forcing)
  use tfm_essentials 
  implicit none

  character(len=*), intent(in)               :: forc_file
  integer, intent(in)                        :: nt
  real(prec), dimension(7,nt), intent(inout) :: forcing

  integer :: n

  open(111, file=forc_file)
    do n = 1, nt, 1
      read(111,*) forcing(2,n), forcing(3,n), forcing(5,n), forcing(4,n)
    end do
  close(111)
end subroutine kohnen03_import_forcing


subroutine kohnen03_sim_info(init_file, forc_file, nz, nt)
  implicit none

  character(len=*), intent(in) :: init_file
  character(len=*), intent(in) :: forc_file

  integer, intent(inout) :: nz
  integer, intent(inout) :: nt

  integer :: stat

  open(111, file=init_file, action='read')
    do
      read(111, fmt='(a)', iostat=stat)
      if ( stat /= 0 ) EXIT
      nz = nz + 1
    end do
  close(111)

  open(222, file=forc_file, action='read')
    do
      read(222, fmt='(a)', iostat=stat)
      if ( stat /= 0 ) EXIT
      nt = nt + 1
    end do
  close(222)
end subroutine kohnen03_sim_info
