program retmip
  use settings
  use tfm_num
  use tfm_tools
  use tfm_output
  use netcdf
  implicit none

  integer, parameter :: NZ_INIT = 600
  integer, parameter :: WRITE_INTER = (24 / 3)

  !integer, parameter :: NT = 1436
  integer, parameter :: NT = 49438
  !integer, parameter :: NT = 1880
  !integer, parameter :: NT = 42901
  !integer, parameter :: NT = 13643

  !character(len=*), parameter :: MODEL = 'Dye-2_16'
  character(len=*), parameter :: MODEL = 'Dye-2_long'
  !character(len=*), parameter :: MODEL = 'FA'
  !character(len=*), parameter :: MODEL = 'Summit'
  !character(len=*), parameter :: MODEL = 'KANU'

  character(len=*), parameter :: INP_FORC = '../RetMIP/RetMIP_forcing_'//MODEL//'.dat'
  character(len=*), parameter :: INP_INIT = '../RetMIP/RetMIP_init_'//MODEL//'.dat'
  character(len=*), parameter :: NC_OUT   = './RetMIP_'//MODEL//'.nc'

  real(prec), dimension(6,NT)             :: forcing
  real(prec), dimension(NT)               :: runoff
  real(prec), dimension(:,:), allocatable :: props

  ! netcdf
  type(out_ncid) :: ncid

  type(sim_models) :: models
  real(prec)       :: dt
  integer          :: np, nz = NZ_INIT
  integer          :: t, n_out = 0
  integer          :: tic, toc
  real(prec)       :: rate

  ! forcing
  call retmip_forcing(NT, INP_FORC, forcing)
  dt = forcing(2,2) - forcing(2,1)

  ! profile initialization
  np = nz + NT
  allocate(props(6,np))
  props = -9999.9
  call retmip_init(NZ, INP_INIT, props(:,1:nz))

  ! model initialization
  call tfm_num_modelinit(                    &
  &  solve_density='arthern2010',            &
  &  solve_temperature='true',               &
  &  solve_heat_capacity='paterson1994',     &
  &  solve_thermal_conductivity='sturm2007', &
  &  solve_liquid='bucket',                  &
  &  models=models                           &
  )

  ! output file initialization
  call tfm_ncout_init(NC_OUT, 0, np, ncid)
  call tfm_ncout_write(ncid, n_out, nz, props, 0.0_prec)

  ! time loop
  call system_clock(tic, rate)
  do t = 1, NT, 1
    call tfm_tools_indicate_tstep(NT, t)

    call tfm_num_surface(       &
    &  nz, np, dt,              &
    &  forcing(:,t),            &
    &  models,                  &
    &  props                    &
    )

    call tfm_num_step(          &
    &  nz, dt,                  &
    &  models=models,           &
    &  props=props(:,1:nz),     &
    &  runoff=runoff(t),        &
    &  liquid_acc=forcing(6,t), &
    &  solid_acc=forcing(5,t)   &
    )
    
    if ( modulo((t - 1), WRITE_INTER) == 0 ) then
      call tfm_ncout_write(ncid, n_out, nz, props, forcing(1,t))
    end if
  end do
  call system_clock(toc, rate)
  write(*, '(a,f10.2,a)') 'time elapsed: ', real(toc - tic) / real(rate), ' s'

  ! close output file
  call tfm_ncout_close(ncid)

  open(111, file=MODEL//'.out', action='write')
  do t = 1, nz, 1
    write(111,*) props(1,t), props(2,t), props(3,t)
  end do
  close(111)

  deallocate(props)
end program retmip


subroutine retmip_init(nz, inp, props)
  use settings
  implicit none

  integer, intent(in)          :: nz
  character(len=*), intent(in) :: inp

  real(prec), dimension(6,nz), intent(inout) :: props

  integer :: n

  open(111, file=inp, action='read')
  
  do n = nz, 1, -1
    read(111,*) props(1,n), props(2,n), props(3,n)
  end do

  close(111)
end subroutine retmip_init


subroutine retmip_forcing(nt, inp, forcing)
  use settings
  use tfm_tools
  implicit none

  integer, intent(in)          :: nt
  character(len=*), intent(in) :: inp

  real(prec), dimension(6,nt), intent(inout) :: forcing

  real(prec), dimension(nt) :: time
  real(prec), dimension(nt) :: melt
  real(prec), dimension(nt) :: acc
  real(prec), dimension(nt) :: temp

  integer    :: n
  real(prec) :: year, month, day

  ! forcing data
  open(111, file=inp, action='read')
  do n = 1, nt, 1
    read(111,*) year, month, day, melt(n), acc(n), temp(n)
    call julian_day(year, month, day, time(n))
  end do
  close(111)

  forcing(1,:) = time                               ! julian day
  forcing(2,:) = (time - time(1)) * (3600.0 * 24.0) ! model time
  forcing(3,:) = temp                               ! temperature
  forcing(4,:) = 315.0                              ! surface density
  forcing(5,:) = acc                                ! solid accumulation
  forcing(6,:) = melt                               ! liquid accumulation
end subroutine retmip_forcing
