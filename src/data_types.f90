module data_types

  use iso_fortran_env, only : dp => real64
  implicit none

  type SU2
     complex(dp), dimension(2,2) :: matrix
  end type SU2

  type(SU2) :: one

contains

  subroutine set_one
    one%matrix = reshape([1.0_dp,0.0_dp,0.0_dp,1.0_dp],[2,2])
  end subroutine set_one
  
end module data_types
