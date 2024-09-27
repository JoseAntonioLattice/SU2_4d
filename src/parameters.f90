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
  logical :: readbeta
  real(dp),  dimension(2) :: b_array 
  integer(i4) :: N_beta
  character(40):: algorithm, smoothing_method

  integer(i4) :: out_smooth_history, master_file_unit

  namelist /input_parameters/ L,Lt,N_thermalization,N_measurements,&
       N_skip,n_time, dt,&
       readbeta,b_array,N_beta,&
       algorithm, smoothing_method

  namelist /lattice_parameters/      L, Lt
  namelist /measurements_parameters/ N_measurements, N_thermalization, N_skip
  namelist /algorithm_parameters/    algorithm
  namelist /beta_parameters/         readbeta,b_array,n_beta
  namelist /smooth_parameters/       smoothing_method,n_time,dt 
  
contains
  
  subroutine read_input()
    !use iso_fortran_env, only : out_unit
    character(256) :: filename
    integer(i4) :: unit
    logical :: lua_condition, smooth_condition

    write(*,'("Please enter the parameters file.")', advance = 'no')
    read(*,'(a)') filename
    write(*,'(" User  typed : ", a)') trim(filename)
    open(newunit = unit, file = trim(filename), status = 'old')

    read(unit, nml = lattice_parameters)
    !write(*,nml = lattice_parameters)

    read(unit, nml = measurements_parameters)
    !write(*,nml = measurements_parameters)

    read(unit, nml = algorithm_parameters)
    !write(*,nml = algorithm_parameters)
    
    read(unit, nml = beta_parameters)
    !write(*,nml = beta_parameters)

    read(unit, nml = smooth_parameters)
    !write(*,nml = smooth_parameters)
    
    close(unit)

    write(*,nml = input_parameters)
    
    lua_condition    = algorithm /= 'metropolis' .or. algorithm /= 'heatbath'
    
    smooth_condition = smoothing_method /= 'cooling'       .or. &
                       smoothing_method /= 'gradient_flow' .or. &
                       smoothing_method /= 'ape_smearing'
    
    if ( L  < 1 )                 stop "L must be > 0"
    if ( Lt < 1 )                 stop "Lt must be > 0"
    if ( N_thermalization < 1 )   stop "N_thermalization must be > 0"
    if ( N_measurements   < 1 )   stop "N_measurements must be > 0"
    if ( N_skip < 1 )             stop "N_skip must be > 0"
    if ( N_time < 1 )             stop "N_time must be > 0"
    if ( dt <= 0.0_dp )           stop "dt must be > 0"
    if ( .not. lua_condition    ) stop "algorithm must be 'metropolis' or 'heatbath'"
    if ( .not. smooth_condition ) stop "smoothing_method must be 'cooling', 'gradient_flow' or 'ape_smearing'"

  end subroutine read_input
end module parameters
