module parameters

  use precision
  implicit none

  integer(i4), parameter :: d = 4
  integer(i4), parameter :: N = 2
  integer(i4) :: L, Lt
  integer(i4) :: N_thermalization
  integer(i4) :: N_measurements
  integer(i4) :: N_skip
  real(dp)    :: bi, bf
  integer(i4) :: N_beta
  character(99):: algorithm
  namelist /input_parameters/ L,Lt,N_thermalization,N_measurements,N_skip,bi,bf,N_beta,algorithm

contains
  
  subroutine read_input()
    character(256) :: filename
    integer(i4) :: unit

    read(*,'(a)') filename
    open(newunit = unit, file = trim(filename), status = 'old')
    read(unit, nml = input_parameters)
    write(*,nml = input_parameters)
  end subroutine read_input
end module parameters
