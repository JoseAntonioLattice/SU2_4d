module parameters

  use precision
  implicit none

  integer(i4), parameter :: d = 4
  integer(i4), parameter :: N = 2
  integer(i4) :: L, Lt
  integer(i4) :: N_thermalization
  integer(i4) :: N_measurements
  integer(i4) :: N_skip
  integer(i4) :: N_time
  real(dp) :: dt
  real(dp),  dimension(2) :: b_array 
  integer(i4) :: N_beta
  character(40):: algorithm, smoothing_method

  integer(i4) :: out_smooth_history

  
  namelist /input_parameters/ L,Lt,N_thermalization,N_measurements,N_skip,n_time, dt,b_array,N_beta,algorithm, smoothing_method

contains
  
  subroutine read_input()
    !use iso_fortran_env, only : out_unit
    character(256) :: filename
    integer(i4) :: unit

    write(*,'("Please enter the parameters file.")', advance = 'no')
    read(*,'(a)') filename
    write(*,'(" User  typed : ", a)') trim(filename)
    open(newunit = unit, file = trim(filename), status = 'old')
    read(unit, nml = input_parameters)
    close(unit)

    if( L < 1 )  stop " L must be > 0"
    if( Lt < 1 ) stop " Lt must be > 0"
    if( N_thermalization < 1 ) stop " N_thermalization  must be > 0"
    if( N_measurements < 1 ) stop " N_measurements must be > 0"
    if( N_skip < 1 ) stop " N_skip must be > 0"
    if( N_time < 1 ) stop " N_time must be > 0"
    if( dt <= 0.0_dp  ) stop " dt must be > 0"
    if( .not. (algorithm /= 'metropolis' .or. algorithm /= 'heatbath')) stop "algorithm must be 'metropolis' or 'heatbath'"
    if( .not. (smoothing_method /= 'cooling' .or. smoothing_method /= 'gradient_flow')) &
         stop "smoothing_method must be 'cooling' or 'gradient_flow'"
    write(*,nml = input_parameters)

  end subroutine read_input
end module parameters
