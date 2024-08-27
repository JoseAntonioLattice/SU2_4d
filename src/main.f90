program SU2_4d
  use parameters
  use arrays
  use dynamics
  use datatypes
  use pbc
  use statistics
  use create_files
  implicit none
  integer :: i, bins
  call read_input()
  allocate(U(d,L,L,L,Lt))
  allocate(beta(n_beta))
  allocate(P(n_measurements))
  call set_periodic_bounds(L,Lt)
  call hot_start(U)
  
  open(unit = 10, file = 'data/Lx='//trim(int2str(L))//'_Lt='&
                         //trim(int2str(Lt))//'_'//trim(algorithm)//'.dat'&
                         ,status = 'unknown')
  do i = 1, N_beta
     beta(i) = bi + (bf - bi)/(N_beta - 1) * (i-1)
     call create_measurements_file(L,Lt,beta(i),algorithm,.true.)
     call thermalization(U,beta(i))
     call measurements(U,beta(i),P)
     call max_jackknife_error_2(P,avr_P,err_P,bins)
     print*, beta(i),avr_P,err_P
     write(10,*) beta(i), avr_P,err_P
     flush(10)
  end do
end program SU2_4d
