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
  call read_beta()
  !allocate(U(d,L,L,L,Lt))
  !print*, 'ok 1'
  !allocate(beta(n_beta))
  !print*, 'ok 2'
  !allocate(P(n_measurements))
  !allocate(Q_den(n_measurements))
  !call create_levicivita()
  !call set_periodic_bounds(L,Lt)
  !call hot_start(U)
  !beta = 0.1_dp
  !open(unit = 10, file = 'data/Lx='//trim(int2str(L))//'_Lt='&
  !                       //trim(int2str(Lt))//'_'//trim(algorithm)//'.dat'&
  !                       ,status = 'unknown')

  !do i = 1, N_beta
  !   !beta(i) = bi + (bf - bi)/(N_beta - 1) * (i-1)
  !   call create_measurements_file(L,Lt,beta(i),algorithm,.true.)
  !   call thermalization(U,beta(i))
  !   call measurements(U,beta(i),P,Q_den)
  !   call max_jackknife_error_2(P,avr_P,err_P,bins)
  !   call max_jackknife_error_2(Q_den%re,avr_Qden%re,err_Qden%re,bins)
  !   call max_jackknife_error_2(Q_den%im,avr_Qden%im,err_Qden%im,bins)
  !   print*, beta(i),avr_P,err_P, avr_Qden%re, err_Qden%re,avr_Qden%im, err_Qden%im
  !   write(10,*) beta(i), avr_P,err_P
  !   flush(10)
  !end do
end program SU2_4d
