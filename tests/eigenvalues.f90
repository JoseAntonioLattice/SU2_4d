program eigenvalues

  use iso_fortran_env, only : dp => real64, i4 => int32
  implicit none

  integer, parameter :: n = 3
  integer, parameter :: lda = 3
  integer, parameter :: ldvl  = n
  integer, parameter :: ldvr = n
  integer, parameter :: lwork = 2*n
  
  complex(dp), dimension(n,n) :: A, W, B, C, expA, one
  complex(dp), dimension(lwork) :: work
  complex(dp), dimension(n) :: eigenv
  complex(dp), dimension(ldvl,n) :: vl
  complex(dp), dimension(ldvr,n) :: vr
  real(dp), dimension(2*n) :: rwork
  integer :: info
  complex(dp) :: det
  W = reshape([0.0_dp,-2.0_dp,0.0_dp,1.0_dp,-3.0_dp,0.0_dp,0.0_dp,0.0_dp,1.0_dp],[n,n])
  one = reshape([1.0_dp,0.0_dp,0.0_dp, 0.0_dp,1.0_dp,0.0_dp, 0.0_dp,0.0_dp,1.0_dp],[n,n])
  W = TA(W,n)
  A = W
  call zgeev('N', 'V', n, A, lda,eigenv, vl, ldvl, vr, ldvr, WORK, lwork, rwork,INFO)
  
  C = vr
  B = 0.0_dp
  B(1,1) = exp(eigenv(1))
  B(2,2) = exp(eigenv(2))
  B(3,3) = exp(eigenv(3))
  det = eigenv(1)*eigenv(2)*eigenv(3)
  !print*, eigenv
  !print*, C
  print*, matmul(matmul(C,B),inv(C,n))
  print*, su3exp(W,n,det)
  !expA = 1/(eigenv(1) - eigenv(2)) * ( exp(eigenv(1)) * (W - eigenv(2) * one ) &
  !      - exp(eigenv(2)) * (W - eigenv(1) * one ))
  !print*, real(expA)
  !print*, 2*exp(-1.0_dp) - exp(-2.0_dp), -2*exp(-1.0_dp) + 2*exp(-2.0_dp),&
  !          exp(-1.0_dp) - exp(-2.0_dp), -exp(-1.0_dp) + 2*exp(-2.0_dp)

contains
  ! -- Returns the inverse of a general squared matrix A
  function inv(A,n) result(Ainv)
    integer(i4), intent(in) :: n
    complex(dp), dimension(n,n), intent(in) :: A
    complex(dp), dimension(n,n) :: Ainv
    complex(dp)            :: work(n)     ! work array for LAPACK                                                                      
    integer         :: info,ipiv(n)     ! pivot indices
    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    !n = 3
    ! SGETRF computes an LU factorization of a general M-by-N matrix A
    
    ! using partial pivoting with row interchanges.
    
    call zGETRF(n,n,Ainv,n,ipiv,info)
    if (info.ne.0) stop 'Matrix is numerically singular!'
    ! SGETRI computes the inverse of a matrix using the LU factorization
    
    ! computed by SGETRF.                                                                                                                    
    call zGETRI(n,Ainv,n,ipiv,work,n,info)
    if (info.ne.0) stop 'Matrix inversion failed!'
  end function inv

  function su3exp(X,n,det)
    integer(i4), intent(in) :: n
    complex(dp), dimension(n,n), intent(in) :: X    
    complex(dp), intent(in) :: det
    complex(dp), dimension(n,n) :: su3exp
    complex(dp) :: t, d,q0,q1,q2,q0_old,q1_old,q2_old
    complex(dp), dimension(n,n) :: X2
    integer(i4) :: Ni, i

    Ni = 14
    
    X2 = matmul(X,X)
    t = tr(X2)/2
    d = det!X(1,1)*X(2,2) - X(1,2)*X(2,1)

    q0_old = 1.0_dp/fact(Ni)
    q1_old = 0.0_dp
    q2_old = q1_old
    
    do i = Ni - 1, 0, -1
       q0 = 1.0_dp/fact(i) + d*q2_old
       q1 = q0_old + t*q2_old
       q2 = q1_old

       q0_old = q0
       q1_old = q1
       q2_old = q2
    end do

    su3exp = q0*one + q1*X + q2*X2
    
  end function su3exp

  function fact(n)
    integer(i4) :: fact, i
    integer(i4), intent(in) :: n

    fact = 1
    if ( n == 0 ) return

    do i = 1, n
       fact = fact * i
    end do
    
  end function fact

  function TA(A,n)
    integer(i4), intent(in) :: n
    complex(dp), dimension(n,n), intent(in) :: A
    complex(dp), dimension(n,n) :: TA

    TA = (A - transpose(conjg(A)))/2 - tr(A - transpose(conjg(A)))/6 * one
  end function TA

  function tr(A)
   complex(dp), dimension(3,3), intent(in) :: A
   complex(dp) :: tr
   tr = A(1,1) + A(2,2) + A(3,3)
  end function tr
  
end program eigenvalues
