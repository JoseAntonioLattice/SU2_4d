program SU2_4d
  use parameters
  use arrays
  use dynamics
  use datatypes
  use pbc
  use statistics
  implicit none
  integer :: i
  call read_input()
  allocate(U(d,L,L,L,Lt))
  allocate(beta(n_beta))
  allocate(P(n_measurements))
  call set_periodic_bounds(L,Lt)
  call hot_start(U)
  print*, det(small_SU2_ran())

  print*, small_SU2_ran()
  open(unit = 10, file = 'data/data.dat', status = 'unknown')
  do i = 1, N_beta
     beta(i) = bi + (bf - bi)/(N_beta - 1) * (i-1)
     call thermalization(U,beta(i))
     call measurements(U,beta(i),P)
     call std_err(P,avr_P,err_P)
     print*, beta(i),avr_P,err_P
     write(10,*) beta(i), avr_P,err_P
     flush(10)
  end do
end program SU2_4d
