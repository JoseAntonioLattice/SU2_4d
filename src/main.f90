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
  use configurations
  implicit none

  integer(i4) ::i ,i_b
  
  call read_input()
  if (readbeta .eqv. .true.)then
     call read_beta()
  else
     call create_beta_array()
  end if
  
  call reserve_memory()
  call create_levicivita()
  call create_one()
  call set_periodic_bounds(L,Lt)

  
  !open(newunit = master_file_unit, file = 'data/Lx='//trim(int2str(L))//'_Lt='&
  !                      //trim(int2str(Lt))//'_'//trim(algorithm)//'.dat'&
  !                       ,status = 'unknown')
  !call simulate_equilibrium(U,beta,master_file_unit,P,Poly,correlation_polyakov_loop)
  !call read_configuration(U,beta(1))

  open(newunit = out_smooth_history, file = 'data/Lx='//trim(int2str(L))//'_Lt='&
                        //trim(int2str(Lt))//'_'//trim(smoothing_method)//'.dat'&
                         ,status = 'unknown')
  
  
  !call hot_start(U)
  !call thermalization(U,beta(1))
  !do i_b = 1, size(beta)
     do i = 1, N_measurements
        !call hot_start(U)
        call read_configuration(U,0.0_dp,i)
        call smooth_configuration(U,0.0_dp,n_time,smoothing_method,out_smooth_history)
     end do
  !end do
  deallocate(U,P,Q_den,Eden,beta)
  
end program SU2_4d
