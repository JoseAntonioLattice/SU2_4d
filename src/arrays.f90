module arrays
  use precision
  use datatypes
  implicit none

  type(SU2), allocatable, dimension(:,:,:,:,:) :: U
  real(dp), allocatable, dimension(:) :: beta, P
  real(dp) :: avr_P, err_P
end module arrays
