module dynamics

  use data_types
  use iso_fortran_env, only : dp => real64, i4 => int32
  implicit none

contains

  subroutine set_memory(U,beta,L,Lt,N_beta,bi,bf)
    type(SU2), allocatable, dimension(:,:,:,:,:) :: U
    real(dp), allocatable, dimension(:) :: beta
    integer(i4), intent(in) :: L,Lt,N_beta
    real(dp), intent(in) :: bi, bf

    integer(i4) :: i
    
    allocate(U(4,L,L,L,Lt))
    allocate(beta(N_beta))

    do i=1,N_beta
       beta(i) = bi + (bf-bi)/(N_beta-1) * (i-1)
    end do
    
  end subroutine set_memory

  subroutine cold_start(U)
    type(SU2), dimension(:,:,:,:,:), intent(out) :: U

    U = one
    
  end subroutine cold_start

  subroutine hot_start(U)
    type(SU2), dimension(:,:,:,:,:), intent(inout) :: U
    integer(i4) :: mu,x,y,z,t
    integer(i4) :: L, Lt

    L = size((U(1,:,1,1,1)))
    Lt = size((U(1,1,1,1,:)))
    do x = 1, L
       do y = 1, L
          do z = 1, L
             do t = 1, Lt
                do mu = 1, 4
                   U(mu,x,y,z,t) = SU2_ran()
                end do
             end do
          end do
       end do
    end do
  end subroutine hot_start

  function SU2_ran()
    type(SU2) :: SU2_ran

    real(dp) :: r(4)
    real(dp), parameter :: eps = 0.5_dp
    complex(dp) :: a,b
    
    call random_number(r)

    r = r - 0.5_dp
    r(1:3) = eps * r(1:3)/norm2(r(1:3))
    r(4) = sign(r(4),r(4)) * sqrt(1.0_dp - eps**2)

    a = cmplx(r(4),r(1),dp)
    b = cmplx(r(3),r(2),dp)

    SU2_ran%matrix = reshape([a,-conjg(b),b,conjg(a)],[2,2]) 
    
  end function SU2_ran
  
end module dynamics
