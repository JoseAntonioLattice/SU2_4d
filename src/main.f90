program SU2_4d
  use parameters
  use arrays
  use dynamics
  use datatypes
  use pbc
  use statistics
  use create_files
  implicit none
  integer :: i, bins, i_t, i_b
  call read_input()
  call read_beta()
  allocate(U(d,L,L,L,Lt))
 
  !allocate(beta(n_beta))
  
  allocate(P(n_measurements,0:n_time))
  allocate(Q_den(n_measurements,0:n_time))
  allocate(Eden(n_measurements,0:n_time))
  call create_levicivita()
  call create_one()
 
  call set_periodic_bounds(L,Lt)
  call hot_start(U)
 
  open(unit = 10, file = 'data/Lx='//trim(int2str(L))//'_Lt='&
                         //trim(int2str(Lt))//'_'//trim(algorithm)//'.dat'&
                         ,status = 'unknown')

  open(unit = 200, file = 'data/wilson.dat', status = 'unknown')
  do i_b =1, size(beta)
     do i = 1, N_measurements
        call hot_start(U)
        call thermalization(U,beta(i_b))
        Eden(i,0) = energy_density(U)
        P(i,0) =  plaquette_value(U)
        q_den(i,0) =  topological_charge(U)
        do i_t = 1, N_time
           call sweeps(U,beta(i_b),'wilson_flow')
           Eden(i,i_t) = energy_density(U)
           P(i,i_t) =  plaquette_value(U)
           q_den(i,i_t) =  topological_charge(U)
        end do
     end do
     do i_t = 0, n_time
        call max_jackknife_error_2(Eden(:,i_t),avr_eden,err_eden,bins)
        call max_jackknife_error_2(P(:,i_t),avr_P,err_P,bins)
        call max_jackknife_error_2(Q_den(:,i_t),avr_Qden,err_Qden,bins)
        write(200,*) i_t*dt,avr_P,err_P,avr_qden,err_Qden,avr_eden,err_eden
     end do
  end do
end program SU2_4d
