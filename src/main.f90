program SU2_4d
  use parameters
  use data_types
  use arrays
  use dynamics
  use statistics
  implicit none

  integer :: i_b
  
  call read_input
  call set_one
  call set_memory(U,beta,L,L,N_beta,bi,bf,plq_action,n_measurements)
  call hot_start(U)

  open( unit = 10, file = 'data/data.dat', status = 'unknown')
  do i_b = 1, n_beta
     call initialization(u,plq_action,beta(i_b),N_thermalization,N_measurements, N_skip)
     call max_jackknife_error_2(plq_action,avr_action,err_action,bins)
     print*,beta(i_b), avr_action, err_action
     write(10,*) beta(i_b), avr_action, err_action
     flush(10)
  end do


  
end program SU2_4d
