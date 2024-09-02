program SU2_4d
  use parameters
  use arrays
  use dynamics
  use datatypes
  use pbc
  use statistics
  use create_files
  implicit none
  integer :: i, bins, i_t
  call read_input()
  call read_beta()
  allocate(U(d,L,L,L,Lt))
  
  !allocate(beta(n_beta))
  
  allocate(P(n_measurements))
  allocate(Q_den(n_measurements))
   allocate(Eden(n_measurements,100))
  call create_levicivita()
  call set_periodic_bounds(L,Lt)
  call hot_start(U)
  !beta = 0.1_dp
  open(unit = 10, file = 'data/Lx='//trim(int2str(L))//'_Lt='&
                         //trim(int2str(Lt))//'_'//trim(algorithm)//'.dat'&
                         ,status = 'unknown')

  open(unit = 200, file = 'data/wilson.dat', status = 'unknown')
  do i = 1, size(beta)
  !   !beta(i) = bi + (bf - bi)/(N_beta - 1) * (i-1)
     call create_measurements_file(L,Lt,beta(i),algorithm,.true.)
     call thermalization(U,beta(i))
     call measurements(U,beta(i),P,Q_den,Eden)
     call max_jackknife_error_2(P,avr_P,err_P,bins)
     call max_jackknife_error_2(Q_den,avr_Qden,err_Qden,bins)
     do i_t = 1, 100
        call max_jackknife_error_2(Eden(:,i_t),avr_eden,err_eden,bins)
        write(200,*) i_t*0.1_dp,avr_eden,err_eden
     end do
     print*, beta(i),avr_P,err_P, avr_Qden, err_Qden
     write(10,*) beta(i), avr_P,err_P
     flush(10)
  end do
end program SU2_4d
