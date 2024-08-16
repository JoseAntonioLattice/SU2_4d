module arrays

  use iso_fortran_env, only : dp => real64, i4 => int32
  use data_types
  implicit none
  
  type(SU2), allocatable, dimension(:,:,:,:,:) :: U
  real(dp), allocatable, dimension(:) :: beta
  real(dp), allocatable, dimension(:) :: plq_action

  real(dp) :: avr_action,err_action
  integer(i4) :: bins

end module arrays
