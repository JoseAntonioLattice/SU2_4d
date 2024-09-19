program eigenvalues

  use iso_fortran_env, only : dp => real64
  implicit none

  integer, parameter :: n = 2
  integer, parameter :: lda = 2
  integer, parameter :: ldvl  = n
  integer, parameter :: ldvr = n
  integer, parameter :: lwork = 2*n
  
  complex(dp), dimension(n,n) :: A, W, B, C
  complex(dp), dimension(lwork) :: work
  complex(dp), dimension(n) :: eigenv
  complex(dp), dimension(ldvl,n) :: vl
  complex(dp), dimension(ldvr,n) :: vr
  real(dp), dimension(2*n) :: rwork
  integer :: info

  W = reshape([1.0_dp,1.0_dp,-1.0_dp,1.0_dp],[2,2])
  
  A = W
  call zgeev('N', 'V', n, A, lda,eigenv, vl, ldvl, vr, ldvr, WORK, lwork, rwork,INFO)
  
  C = vr
  B = 0.0_dp
  B(1,1) = eigenv(1)
  B(2,2) = eigenv(2)

  print*, eigenv
  print*, C
  print*, matmul(matmul(C,B),inv(C))

contains
  ! -- Returns the inverse of a general squared matrix A
  function inv(A) result(Ainv)
    complex(dp), dimension(n,n), intent(in) :: A
    complex(dp), dimension(n,n) :: Ainv
    complex(dp)            :: work(2)     ! work array for LAPACK                                                                      
    integer         :: n,info,ipiv(2)     ! pivot indices
    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    n = 2
    ! SGETRF computes an LU factorization of a general M-by-N matrix A
    
    ! using partial pivoting with row interchanges.
    
    call zGETRF(n,n,Ainv,n,ipiv,info)
    if (info.ne.0) stop 'Matrix is numerically singular!'
    ! SGETRI computes the inverse of a matrix using the LU factorization
    
    ! computed by SGETRF.                                                                                                                    
    call zGETRI(n,Ainv,n,ipiv,work,n,info)
    if (info.ne.0) stop 'Matrix inversion failed!'
  end function inv
  
end program eigenvalues
