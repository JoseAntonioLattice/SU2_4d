module arrays
  use precision
  use datatypes
  implicit none

  type(SU2), allocatable, dimension(:,:,:,:,:) :: U
  real(dp), allocatable, dimension(:) :: beta, P
  complex(dp),allocatable, dimension(:) :: Q_den 
  real(dp) :: avr_P, err_P
  complex(dp) :: avr_Qden, err_Qden
end module arrays
