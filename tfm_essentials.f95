module tfm_essentials
  use tfm_constants
  implicit none

  contains

  subroutine tfm_essentials_do_nothing(nz, variable)
    implicit none

    integer, intent(in)                   :: nz
    real(prec), dimension(nz), intent(in) :: variable
    real(prec), dimension(nz)             :: nothing
!-----------------------------------------------------------------------
! Subroutne: tfm_essentials_do_nothing
!
! Routine to do "nothing" with a given variable. Used to avoid warning
! at compile time.
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of the given variable.
!   variable: Variable to do nothing about.
!-----------------------------------------------------------------------
    
    nothing = variable
  end subroutine tfm_essentials_do_nothing
  

  subroutine tfm_essentials_mean_acc(nz, depth, density, age, mean_acc)
    use tfm_constants
    implicit none

    integer, intent(in)                      :: nz
    real(prec), dimension(nz), intent(in)    :: depth
    real(prec), dimension(nz), intent(in)    :: density
    real(prec), dimension(nz), intent(in)    :: age
    real(prec), dimension(nz), intent(inout) :: mean_acc

    integer :: n

!-----------------------------------------------------------------------
! Subroutine : tfm_essentials_mean_acc
!
! The subroutine computes the mean accumulation rate from the age of
! the firn profile. The concept follow the idea of calculating the
! mean accumulation rate over the life time of a firn parcel.
!
! See for example:
! Stevens, C. M., Verjans, V., Luding, J. M. D., Kahle, E. C.,
! Horlings, A. N., Horlings, B. I., and Waddington, E. D. The Community
! Firn Model (CFM) v1.0. Geosci. Model. Dev., 13, 4355-4377, (2020).
! https://doi.org/10.5194/gmd-13-4355-2020
!
! Author: Timm Schultz
!
! Arguments:
!   nz: Dimension of variables "depth", "density", "age", "mean_acc".
!   depth: Depth of the firn profile (m).
!   density: Density of the firn profile (kg m**-3).
!   age: Age of the firn profile (s).
!   mean_acc - on input: Variable to store the mean accumulation rate.
!
! Result:
!   mean_acc - on output: Mean accumulation rate (m weq. a**-1).
!-----------------------------------------------------------------------
    
    mean_acc(nz) = 0.0
    do n = nz - 1, 1, -1
      mean_acc(n) = (depth(n+1) - depth(n)) * (density(n) / WATER_DENSITY)
      mean_acc(n) = mean_acc(n) + mean_acc(n+1)
    end do
  
    mean_acc(1:nz-1) = mean_acc(1:nz-1) / (age(1:nz-1) / SECONDS_YEAR)
  end subroutine tfm_essentials_mean_acc
end module tfm_essentials


module tfm_llStructure
  use tfm_constants
  implicit none


  type llNode
    real(prec), dimension(:), allocatable :: data
    type(llNode), pointer                 :: next
  end type llNode


  type linkedList
    integer :: node_size = 1000
    integer :: hind = 1
    integer :: tind = 1001
    integer :: length = 0

    type(llNode), pointer :: head => null()
    type(llNode), pointer :: tail => null()
  end type linkedList


  type llProps
    type(linkedList) :: depth
    type(linkedList) :: density
    type(linkedList) :: temperature
    type(linkedList) :: heatcap
    type(linkedList) :: thermcond
    type(linkedList) :: liquidwater
    type(linkedList) :: age
    type(linkedList) :: grain_radius
  end type llProps


  contains


  subroutine llSetNodeSize(self, node_size)
    implicit none

    type(linkedList), intent(inout) :: self
    integer, intent(in)             :: node_size

    self%node_size = node_size
    self%tind = (self%tind + 1)
  end subroutine llSetNodeSize


  subroutine llAppendNode(self)
    implicit none

    type(linkedList), intent(inout) :: self
    type(llNode), pointer           :: new_node

    ! generate new node
    allocate(new_node)
    allocate(new_node%data(self%node_size))
    new_node%next => null()

    ! if the appended node is the first node, head and tail are the new node
    if ( .not. associated(self%head) ) then
      self%head => new_node
      self%tail => new_node

    ! otherwise the new node bcomes the tail
    else
      self%tail%next => new_node
      self%tail      => new_node

    end if

    ! change linked list attributes accordingly
    self%tind = 1
  end subroutine llAppendNode


  subroutine llAppendData(self, n, data)
    implicit none

    type(linkedList), intent(inout)      :: self
    integer, intent(in)                  :: n
    real(prec), dimension(n), intent(in) :: data

    integer :: remainder
    integer :: free
    integer :: k, l, i, j

    remainder = n
    i = 1

    do while ( remainder > 0 )
      
      ! if necessary append another node
      if ( self%tind > self%node_size ) then
        call llAppendNode(self)
      end if

      ! node and data indices
      free = (self%node_size - self%tind + 1)
      k = self%tind
      l = min((k + free - 1), (k + remainder - 1))
      j = (i + (l - k))

      ! store data
      self%tail%data(k:l) = data(i:j)

      ! change counting variables accordingly
      self%tind = (self%tind + (l - k) + 1)
      i = (i + (l - k) + 1)
      remainder = (remainder - free)
    end do

    self%length = (self%length + n)
  end subroutine llAppendData


  subroutine llDropHeadNode(self)
    implicit none

    type(linkedList), intent(inout) :: self
    type(llNode), pointer           :: former_head_node

    former_head_node => self%head
    self%head => self%head%next

    deallocate(former_head_node%data)
    nullify(former_head_node)
  end subroutine llDropHeadNode


  subroutine llDropTailNode(self)
    implicit none

    type(linkedList), intent(inout) :: self
    type(llNode), pointer           :: new_tail_node

    new_tail_node => self%head
    do while ( associated(new_tail_node%next%next) )
      new_tail_node => new_tail_node%next
    end do

    deallocate(self%tail%data)
    self%tail => new_tail_node
    nullify(self%tail%next)
  end subroutine llDropTailNode


  subroutine llDropData(self, n)
    implicit none

    type(linkedList), intent(inout) :: self
    integer, intent(in)             :: n

    if ( n > 0 ) then
      self%hind = (self%hind + n)

      ! if the head iis empty delete the current head node
      do while ( self%hind >= self%node_size )
        call llDropHeadNode(self)
        self%hind = (self%hind - self%node_size)
      end do

    else if ( n < 0 ) then
      self%tind = (self%tind + n)

      ! if the tail is empty delete the current tail node
      do while ( self%tind < 1 )
        call llDropTailNode(self)
        self%tind = (self%tind + self%node_size)
      end do

    else
      continue
    end if

    self%length = (self%length - abs(n))
  end subroutine llDropData


  function llGetData(self) result(output)
    implicit none

    type(linkedList), intent(in)       :: self
    real(prec), dimension(self%length) :: output

    type(llNode), pointer :: curr_node
    integer               :: k, l, i, j

    ! data from the head
    i = 1
    j = min((self%node_size - self%hind + 1), self%length)
    k = self%hind
    l = min(self%node_size, (self%hind + self%length - 1))
    output(i:j) = self%head%data(k:l)

    if ( associated(self%head, self%tail) ) RETURN

    ! data from the middle part
    curr_node => self%head%next
    do while ( associated(curr_node%next) )
      i = (j + 1)
      j = (i + self%node_size - 1)
      output(i:j) = curr_node%data(1:l)

      curr_node => curr_node%next
    end do

    ! data from the tail
    i = (j + 1)
    j = (i + (self%tind - 1) - 1)
    l = (self%tind - 1)
    output(i:j) = self%tail%data(1:l)
  end function llGetData


  subroutine llUpdateList(self, data)
    implicit none

    type(linkedList), intent(inout)      :: self
    real(prec), dimension(self%length), intent(in) :: data

    integer               :: k, l, i, j
    type(llNode), pointer :: curr_node

    ! head
    i = 1
    j = min((self%node_size - self%hind + 1), self%length)
    k = self%hind
    l = min(self%node_size, (self%hind + self%length - 1))
    self%head%data(k:l) = data(i:j)

    if ( associated(self%head, self%tail) ) RETURN

    ! middle part
    curr_node => self%head%next
    do while ( associated(curr_node%next) )
      i = (j + 1)
      j = (i + self%node_size - 1)
      curr_node%data(:) = data(i:j)

      curr_node => curr_node%next
    end do

    ! tail
    i = (j + 1)
    j = (i + (self%tind - 1) - 1)
    k = 1
    l = (self%tind - 1)
    self%tail%data(k:l) = data(i:j)
  end subroutine llUpdateList


  function llGetFirst(self) result(first)
    implicit none

    type(linkedList), intent(in) :: self
    real(prec)                   :: first

    first = self%head%data(self%hind)
  end function llGetFirst


  function llGetLast(self) result(last)
    implicit none

    type(linkedList), intent(in)       :: self
    real(prec)                         :: last
    real(prec), dimension(self%length) :: all_data

    if ( (self%tind - 1) /= 0 ) then
      last = self%tail%data((self%tind - 1))

    ! if the last index is the first of the tail
    ! get the last entrace of all data
    else
      all_data = llGetData(self)
      last = all_data(self%length)

    end if
  end function llGetLast


  subroutine llFreeList(self)
    implicit none

    type(linkedList), intent(inout) :: self
    type(llNode), pointer           :: curr_node
    type(llNode), pointer           :: next_node

    curr_node => self%head
    do while ( associated(curr_node) )
      next_node => curr_node%next

      deallocate(curr_node%data)
      nullify(curr_node)

      curr_node => next_node
    end do
  end subroutine llFreeList


  subroutine llPropsFree(self)
    implicit none

    type(llProps), intent(inout) :: self

    call llFreeList(self%depth)
    call llFreeList(self%density)
    call llFreeList(self%grain_radius)
    call llFreeList(self%temperature)
    call llFreeList(self%heatcap)
    call llFreeList(self%thermcond)
    call llFreeList(self%liquidwater)
    call llFreeList(self%age)
  end subroutine llPropsFree

  subroutine llPropsDropData(self, n)
    implicit none

    type(llProps), intent(inout) :: self
    integer, intent(in)          :: n

    call llDropData(self%depth,        n)
    call llDropData(self%density,      n)
    call llDropData(self%grain_radius, n)
    call llDropData(self%temperature,  n)
    call llDropData(self%heatcap,      n)
    call llDropData(self%thermcond,    n)
    call llDropData(self%liquidwater,  n)
    call llDropData(self%age,          n)
  end subroutine llPropsDropData
end module tfm_llStructure
