program binary_check

  implicit none

  integer :: m ,n ,o, p, r
  complex(8), allocatable, dimension(:,:,:,:,:) :: a, b, c

  m = 4
  n = 10
  o = 10
  p = 10
  r = 10

  
  print*, "The size of the storage will be: ", 3*m*n*o*p*r*16, " bytes."
  allocate(a(m,n,o,p,r),b(m,n,o,p,r),c(m,n,o,p,r))

  a = 2.0d0
  b = 10.0d0
  c = a + b

  open(unit = 69, file = 'conf.bin', form = 'unformatted', access = "sequential")
  write(69) a
  write(69) b
  write(69) c
  close(69)
  
end program binary_check
  
