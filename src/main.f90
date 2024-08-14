program SU2_4d
  use parameters
  use data_types
  use arrays
  use dynamics
  implicit none

  call read_input
  call set_one
  call set_memory(U,beta,L,L,N_beta,bi,bf)
  call hot_start(U)
  print*, beta
  print*, U(1,1,1,1,1)
end program SU2_4d
