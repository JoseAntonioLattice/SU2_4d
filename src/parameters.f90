module parameters

  use iso_fortran_env, only : dp => real64, i4 => int32
  implicit none

  integer(i4) :: L
  integer(i4) :: N_measurements
  integer(i4) :: N_skip
  integer(i4) :: N_thermalization
  integer(i4) :: N_beta

  real(dp) :: bi, bf

  namelist /input_parameters/ L,N_measurements,N_skip,N_thermalization,N_beta,bi,bf

contains

  subroutine read_input
    integer(i4) :: inunit
    character(256) :: parameters_file

    write(*,'(a)') 'Please, type the parameters file name:'
    read(*,'(a)') parameters_file
    write(*,'(2a)') 'User typed: ', trim(parameters_file)

    open(newunit = inunit, file = trim(parameters_file), status = 'old')
    read(inunit, nml = input_parameters)
    close(inunit)
    write(*,nml = input_parameters)
    
  end subroutine read_input

end module parameters
