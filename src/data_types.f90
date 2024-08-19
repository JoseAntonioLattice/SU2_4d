module data_types

  use iso_fortran_env, only : dp => real64
  implicit none

  type SU2
     complex(dp), dimension(2,2) :: matrix
  end type SU2

  type(SU2) :: one


  interface operator(+)
     module procedure mat_sum
  end interface operator(+)

  interface operator(-)
     module procedure mat_sub
  end interface operator(-)

  interface operator(*)
     module procedure mat_mul
  end interface operator(*)

  interface operator(/)
     module procedure mat_divscalar
  end interface operator(/)

contains

  pure function mat_sum(a,b) result(c)
    type(SU2), intent(in) :: a, b
    type(SU2) :: c
    c%matrix = a%matrix + b%matrix
  end function mat_sum

  pure function mat_sub(a,b) result(c)
    type(SU2), intent(in) :: a, b
    type(SU2) :: c
    c%matrix = a%matrix - b%matrix
  end function mat_sub

  pure function mat_mul(a,b) result(c)
    type(SU2), intent(in) :: a, b
    type(SU2) :: c
    !c%matrix = matmul(a%matrix,b%matrix)
    c%matrix(1,1) = a%matrix(1,1)*b%matrix(1,1) - a%matrix(1,2)*conjg(b%matrix(1,2))
    c%matrix(1,2) = a%matrix(1,1)*b%matrix(1,2) + a%matrix(1,2)*conjg(b%matrix(1,1))
    c%matrix(2,1) = -conjg(c%matrix(1,2))
    c%matrix(2,2) =  conjg(c%matrix(1,1))

  end function mat_mul

  pure function mat_divscalar(A,s) result(C)
    type(SU2), intent(in) :: A
    real(dp), intent(in) :: s
    type(SU2) :: C

    C%matrix = A%matrix/s
  end function mat_divscalar
  
  subroutine set_one
    one%matrix = reshape([1.0_dp,0.0_dp,0.0_dp,1.0_dp],[2,2])
  end subroutine set_one

  pure function dagger(U) result(U_res)
    type(SU2), intent(in) :: U
    type(SU2) :: U_res

    U_res%matrix = transpose(conjg(U%matrix))
  end function dagger

  pure function tr(U) result(trace)
    type(SU2), intent(in) :: U
    complex(dp) :: trace
    

    trace = U%matrix(1,1) + U%matrix(2,2) 

  end function tr

  function det(U)
    type(SU2), intent(in) :: U
    real(dp) :: det

    det = abs(u%matrix(1,1))**2 + abs(u%matrix(1,2))**2
    
  end function det
  
end module data_types
