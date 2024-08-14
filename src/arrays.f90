module arrays

  use iso_fortran_env, only : dp => real64, i4 => int32
  use data_types
  implicit none
  
  type(SU2), allocatable, dimension(:,:,:,:,:) :: U
  real(dp), allocatable, dimension(:) :: beta

end module arrays
