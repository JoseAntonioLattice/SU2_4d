module local_update_algorithms
  use precision
  use datatypes
  use observables, only : DS
  use wilson_loops
  implicit none
  real(dp), parameter :: pi = acos(-1.0_dp)
  
contains

  subroutine metropolis(U,x,mu,beta)
    type(SU2), dimension(:,:,:,:,:), intent(inout) :: U
    integer(i4), intent(in) :: x(4), mu
    real(dp), intent(in) :: beta
    type(SU2) :: Up
    real(dp) :: r, prob, DeltaS

    Up = SU2_ran()
    DeltaS = DS(U,Up,x,mu,beta)
    prob = min(1.0_dp,exp(-DeltaS))

    call random_number(r)
    if( prob >= r ) U(mu,x(1),x(2),x(3),x(4)) = Up

  end subroutine metropolis

  subroutine heatbath(U,x,mu,beta)
    type(SU2), dimension(:,:,:,:,:), intent(inout) :: U
    integer(i4), intent(in) :: x(4), mu
    real(dp), intent(in) :: beta
    type(SU2) :: sigma,V,Xmat
    real(dp) :: alpha
    logical :: condition

    real(dp) :: r(3), lambdasq, s,n(0:3), W(2)
    complex(dp) :: a,b

    sigma = staples(U,x,mu)
    alpha = sqrt(det(sigma))

    V = sigma/alpha

    condition = .false.
    do while (condition .eqv. .false.)
       call random_number(r)
       r = 1.0_dp - r

       lambdasq = -1/(2*alpha*beta) * (log(r(1)) + (cos(2*pi*r(2)))**2 * log(r(3)))

       call random_number(s)
       if (s**2 <= 1.0_dp - lambdasq) condition = .true.
    end do
    n(0) = 1.0_dp - 2*lambdasq

    call random_number(W)
    n(1) = 1.0_dp - 2*W(1)
    n(2) = sqrt(1.0_dp - (n(1))**2) * cos(2*pi*W(2))
    n(3) = sqrt(1.0_dp - (n(1))**2) * sin(2*pi*W(2))

    n(1:3) = sqrt(1.0_dp - (n(0))**2) * n(1:3)

    a = cmplx(n(0),n(1),dp)
    b = cmplx(n(2),n(3),dp)
    Xmat = SU2_matrix(a,b)

    U(mu,x(1),x(2),x(3),x(4)) = Xmat * V

  end subroutine heatbath

  subroutine overrelaxation(U,x,mu)
    type(SU2), dimension(:,:,:,:,:), intent(inout) :: U
    integer(i4), intent(in) :: x(4), mu
    type(SU2) :: sigma,V
    real(dp) :: alpha
    
    sigma = staples(U,x,mu)
    alpha = sqrt(det(sigma))
    V = sigma/alpha
    U(mu,x(1),x(2),x(3),x(4)) = V * dagger(U(mu,x(1),x(2),x(3),x(4))) * V
    !U(mu,x(1),x(2),x(3),x(4)) = tr(U(mu,x(1),x(2),x(3),x(4))*dagger(V))* V - U(mu,x(1),x(2),x(3),x(4))
  end subroutine overrelaxation
  

  
end module local_update_algorithms
