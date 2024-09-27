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

  
  open(newunit = master_file_unit, file = 'data/Lx='//trim(int2str(L))//'_Lt='&
                        //trim(int2str(Lt))//'_'//trim(algorithm)//'.dat'&
                         ,status = 'unknown')
  call simulate_equilibrium(U,beta,master_file_unit,P,Poly,correlation_polyakov_loop)
  deallocate(U,P,Q_den,Eden,beta)
  
end program SU2_4d
