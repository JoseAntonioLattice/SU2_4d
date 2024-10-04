module SU2_math
  
  use iso_fortran_env, only : dp => real64, i4 => int32
  implicit none
  complex(dp), parameter :: i = (0.0_dp,1.0_dp)
  
  type SU2
     complex(dp), dimension(2,2) :: matrix
  end type SU2

  type su2_algebra
     real(dp), dimension(3) :: elements
  end type su2_algebra


contains

  function SU2_ran()
    type(SU2) :: SU2_ran
    real(dp) :: r(4)
    complex(dp) :: a, b
    
    call random_number(r)
    r = r/norm2(r)

    a = cmplx(r(1),r(2),dp)
    b = cmplx(r(3),r(4),dp)
    
    SU2_ran%matrix(1,1) = a
    SU2_ran%matrix(2,1) = -conjg(b)
    SU2_ran%matrix(1,2) = b
    SU2_ran%matrix(2,2) =  conjg(a)

  end function SU2_ran

  function SU2_alg_ran()
    type(SU2) :: SU2_alg_ran
    real(dp) :: r(3)
    complex(dp) :: z
    real(dp) :: a
    
    call random_number(r)
    r = 100*r

    z = cmplx(r(1),r(2),dp)
    a = r(3)

    print*, r
    SU2_alg_ran%matrix(1,1) =  a
    SU2_alg_ran%matrix(2,1) =  z
    SU2_alg_ran%matrix(1,2) = conjg(z)
    SU2_alg_ran%matrix(2,2) = -a

  end function SU2_alg_ran
  
  function tr(V)
    type(SU2)   :: V 
    real(dp) :: tr

    tr = 2*real(V%matrix(1,1))! + V%matrix(2,2) 
  end function tr

  
end module SU2_math

program liealgebra
  use SU2_math
  implicit none

  type(SU2) :: T1, T2, T3, U, W1, W2, W3
  
  T1%matrix = reshape([(0.0_dp,0.0_dp),(1.0_dp,0.0_dp),(1.0_dp,0.0_dp),(0.0_dp,0.0_dp)], [2,2])/2
  T2%matrix = reshape([(0.0_dp,0.0_dp),i,-i,(0.0_dp,0.0_dp)], [2,2])/2
  T3%matrix = reshape([(1.0_dp,0.0_dp),(0.0_dp,0.0_dp),(0.0_dp,0.0_dp),(-1.0_dp,0.0_dp)], [2,2])/2

  W1%matrix = matmul(T1%matrix,T1%matrix)
  print*, tr(W1)

  W1%matrix = matmul(T2%matrix,T2%matrix)
  print*, tr(W1)

  W1%matrix = matmul(T3%matrix,T3%matrix)
  print*, tr(W1)

  
  U = su2_alg_ran()
  W1%matrix = matmul(U%matrix,T1%matrix)
  W2%matrix = matmul(U%matrix,T2%matrix)
  W3%matrix = matmul(U%matrix,T3%matrix)
  print*, tr(W1), tr(W2),tr(W3)
  

end program liealgebra























































