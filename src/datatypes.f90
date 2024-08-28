module datatypes

  use precision
  use parameters, only : N
  implicit none

  type SU2
     complex(dp), dimension(N,N) :: matrix
  end type SU2

  interface operator(+)
     module procedure :: mat_sum
  end interface operator(+)

  interface operator(-)
     module procedure :: mat_sub
  end interface operator(-)

  interface operator(*)
     module procedure :: mat_mul
  end interface operator(*)

  interface operator(*)
     module procedure :: mat_mul_scalar
  end interface operator(*)

  
  interface operator(/)
     module procedure :: mat_div
  end interface operator(/)
contains

  pure function mat_sum(A,B) result(C)
    type(SU2), intent(in) :: A,B
    type(SU2) :: C
    C%matrix = A%matrix + B%matrix
  end function mat_sum

  pure function mat_sub(A,B) result(C)
    type(SU2), intent(in) :: A,B
    type(SU2) :: C
    C%matrix = A%matrix - B%matrix
  end function mat_sub

  pure function mat_mul_scalar(a,B) result(C)
    real(dp), intent(in) :: a
    type(SU2), intent(in) :: B
    type(SU2) :: C
    C%matrix = A*B%matrix
  end function mat_mul_scalar

  
  pure function mat_mul(a,b) result(c)
    type(SU2), intent(in) :: a, b
    type(SU2) :: c
  
    c%matrix(1,1) = a%matrix(1,1)*b%matrix(1,1) - a%matrix(1,2)*conjg(b%matrix(1,2))
    c%matrix(1,2) = a%matrix(1,1)*b%matrix(1,2) + a%matrix(1,2)*conjg(b%matrix(1,1))
    c%matrix(2,1) = -conjg(c%matrix(1,2))
    c%matrix(2,2) =  conjg(c%matrix(1,1))

  end function mat_mul

  pure function mat_div(A,b) result(C)
    type(SU2), intent(in) :: A
    real(dp), intent(in) :: b
    type(SU2) :: C
    C%matrix = A%matrix/b
  end function mat_div
  
  function det(U)
    type(SU2), intent(in) :: U
    real(dp) :: det
    det = (abs(U%matrix(1,1)))**2 + (abs(U%matrix(1,2)))**2 
  end function det

  pure function tr(U)
    type(SU2), intent(in) :: U
    real(dp) :: tr
    tr = 2*real(U%matrix(1,1))
  end function tr

  pure function dagger(U)
    type(SU2), intent(in) :: U
    type(SU2) :: dagger
    dagger%matrix = transpose(conjg(U%matrix))
  end function dagger
  
end module datatypes
