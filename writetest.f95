program writetest
  implicit none

  print *, 'test'
  write(*, '(a)', advance='no') '\b\b'
end program writetest
