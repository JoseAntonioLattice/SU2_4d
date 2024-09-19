module datatypes

  use precision
  use parameters, only : N
  implicit none

  type SU2
     complex(dp), dimension(N,N) :: matrix
  end type SU2

  type(SU2) :: one
  integer(i4), dimension(4,4,4,4) :: levi_civita
  integer(i4), dimension(24,4) :: levi_civita_indices
  
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

  interface tr
     module procedure :: tr_su2
     module procedure :: tr_re
  end interface tr

contains
  
  subroutine create_one()
    one%matrix = reshape([(1.0_dp,0.0_dp),(0.0_dp,0.0_dp),(0.0_dp,0.0_dp),(1.0_dp,0.0_dp)],[2,2])
  end subroutine create_one

  subroutine create_levicivita
    integer(i4) :: i,j,k,l,ii
    ii = 0
    do i = 1, 4
       do j = 1, 4
          do k = 1, 4
             do l = 1, 4
                levi_civita(i,j,k,l) = feps([i,j,k,l])
                if (levi_civita(i,j,k,l) /= 0)then
                   ii = ii + 1
                   !print*, i,j,k,l,levi_civita(i,j,k,l),ii
                   levi_civita_indices(ii,:) = [i,j,k,l]
                end if
             end do
          end do
       end do
    end do

  end subroutine create_levicivita

  pure function feps(x) result(f)
    integer(i4) :: f
    integer(i4), intent(in) :: x(4)
    integer(i4) :: i1, i2
    f = 1
    do i1 = 1,  3
       do i2 = i1 + 1, 4
          f = f * sgn_int(x(i2) - x(i1))
       end do
    end do
    
  end function feps

  pure function sgn(x)
    real(dp), intent(in) :: x
    integer(i4) :: sgn

    if( x > 0.0_dp )then
       sgn = 1
    elseif( x < 0.0_dp)then
       sgn = -1
    else
       sgn = 0
    end if

  end function sgn

  pure function sgn_int(x)
    integer(i4), intent(in) :: x
    integer(i4) :: sgn_int

    if( x > 0 )then
       sgn_int = 1
    elseif( x < 0)then
       sgn_int = -1
    else
       sgn_int = 0
    end if

  end function sgn_int

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

  pure function tr_su2(U)
    type(SU2), intent(in) :: U
    real(dp) :: tr_su2
    tr_su2 = 2*real(U%matrix(1,1))
  end function tr_su2

  pure function tr_re(U)
    real(dp),dimension(2,2), intent(in) :: U
    real(dp) :: tr_re
    tr_re = U(1,1) + U(2,2)
  end function tr_re
  
  pure function dagger(U)
    type(SU2), intent(in) :: U
    type(SU2) :: dagger
    dagger%matrix = transpose(conjg(U%matrix))
  end function dagger

  function mat_exp(W) result(res)
    type(SU2), intent(in) :: W
    type(SU2) :: res, B, C

    integer, parameter :: n = 2
    integer, parameter :: lda = 2
    integer, parameter :: ldvl  = n
    integer, parameter :: ldvr = n
    integer, parameter :: lwork = 2*n

    complex(dp), dimension(n,n) :: A
    complex(dp), dimension(lwork) :: work
    complex(dp), dimension(n) :: eigenv
    complex(dp), dimension(ldvl,n) :: vl
    complex(dp), dimension(ldvr,n) :: vr
    real(dp), dimension(2*n) :: rwork
    integer :: info

    A = W%matrix
    call zgeev('N', 'V', n, A, lda,eigenv, vl, ldvl, vr, ldvr, WORK, lwork, rwork,INFO)

    C%matrix = vr
    B%matrix = 0.0_dp
    B%matrix(1,1) = exp(eigenv(1))
    B%matrix(2,2) = exp(eigenv(2))
 
    res = C*B*inv(C)

  end function mat_exp

  ! -- Returns the inverse of a general squared matrix A
  function inv(A) result(Ainv)
    implicit none
    type(SU2), intent(in)::  A
    type(SU2) :: Ainv
    complex(dp)            :: work(2)     ! work array for LAPACK                                                                      
    integer         :: n,info,ipiv(2)     ! pivot indices
    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv%matrix = A%matrix
    n = 2
    ! SGETRF computes an LU factorization of a general M-by-N matrix A
    
    ! using partial pivoting with row interchanges.
    
    call zGETRF(n,n,Ainv%matrix,n,ipiv,info)
    if (info.ne.0) stop 'Matrix is numerically singular!'
    ! SGETRI computes the inverse of a matrix using the LU factorization
    
    ! computed by SGETRF.                                                                                                                    
    call zGETRI(n,Ainv%matrix,n,ipiv,work,n,info)
    if (info.ne.0) stop 'Matrix inversion failed!'
  end function inv
  
  function SU2_ran() 
    type(SU2) :: SU2_ran
    real(dp) :: r(4)
    complex(dp) :: a,b

    call random_number(r)
    r = r - 0.5_dp
    r = r / norm2(r)
    a = cmplx(r(1),r(2),dp)
    b = cmplx(r(3),r(4),dp)
    SU2_ran = SU2_matrix(a,b)
    
  end function SU2_ran

  function small_SU2_ran() 
    type(SU2) :: small_SU2_ran
    real(dp) :: r(0:3)
    complex(dp) :: a,b
    real(dp), parameter :: eps = 0.1_dp

    call random_number(r)
    r = r - 0.5_dp
    r(1:3) = eps * r(1:3) / norm2(r(1:3))
    r(0) = sgn(r(0)) * sqrt(1.0_dp - eps**2)
    a = cmplx(r(0),r(1),dp)
    b = cmplx(r(2),r(3),dp)
    small_SU2_ran = SU2_matrix(a,b)
    
  end function small_SU2_ran

  pure function SU2_matrix(a,b)
    type(SU2) :: SU2_matrix
    complex(dp), intent(in) :: a,b

    SU2_matrix%matrix(1,1) = a
    SU2_matrix%matrix(1,2) = b
    SU2_matrix%matrix(2,1) = -conjg(b)
    SU2_matrix%matrix(2,2) =  conjg(a)
  end function SU2_matrix

end module datatypes
