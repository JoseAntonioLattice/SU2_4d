program SU2_4d

  use parameters
  use arrays
  use dynamics
  use datatypes
  use pbc
  use statistics
  use create_files
  use observables
  use starts
  use number2string_mod

  implicit none

  integer :: i, bins, i_t, i_b
  real(dp) :: t2,t1


  call read_input()
  call read_beta()
  call reserve_memory()
  call create_levicivita()
  call create_one()
  call set_periodic_bounds(L,Lt)

  !go to 200
  
  open(unit = 10, file = 'data/Lx='//trim(int2str(L))//'_Lt='&
                         //trim(int2str(Lt))//'_'//trim(algorithm)//'.dat'&
                         ,status = 'unknown')

  open(unit = 200, file = 'data/'//trim(int2str(L))//trim(smoothing_method)//'.dat', status = 'unknown')
  open(newunit = out_smooth_history, file = 'data/'//trim(int2str(L))// &
       trim(smoothing_method)//'_history.dat', status = 'unknown')

  call cpu_time(t1)
  call progress_bar(0.0)
  do i_b =1, size(beta)
     do i = 1, N_measurements
        call hot_start(U)
        call thermalization(U,beta(i_b))
        Eden(i,0) = E(U,'plaquette')
        P(i,0) =  plaquette_value(U)
        q_den(i,0) =  topological_charge(U,'plaquette')
        call smooth_configuration(U,beta(i_b),n_time,eden(i,:),P(i,:),q_den(i,:),smoothing_method,out_smooth_history)
        call progress_bar(i/real(N_measurements))
     end do
     do i_t = 0, n_time
        call max_jackknife_error_2(Eden(:,i_t),avr_eden,err_eden,bins)
        call max_jackknife_error_2(P(:,i_t),avr_P,err_P,bins)
        call max_jackknife_error_2(Q_den(:,i_t),avr_Qden,err_Qden,bins)
        write(200,*) i_t,avr_P,err_P,avr_qden,err_Qden,avr_eden,err_eden
     end do
  end do
  deallocate(U,P,Q_den,Eden)
  call cpu_time(t2)
  write(out_smooth_history,*) "Time : ", t2-t1 , "secs" 
  200 print*, "Fin."
end program SU2_4d
