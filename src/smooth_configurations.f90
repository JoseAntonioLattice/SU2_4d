module smooth_configurations
  use precision
  use datatypes
  use wilson_loops
  implicit none

contains

  subroutine cooling(U,x,mu)
    type(SU2), dimension(:,:,:,:,:), intent(inout) :: U
    integer(i4), intent(in) :: mu, x(4)
    type(SU2) :: sigma
    real(dp) :: alpha
    sigma = staples(U,x,mu)
    alpha = sqrt(det(sigma))
    U(mu,x(1),x(2),x(3),x(4)) = sigma/alpha
  end subroutine cooling
  
  subroutine gradient_flow(U,x,mu)
    use parameters, only : dt
    integer(i4), intent(in) :: x(4), mu
    type(SU2), dimension(:,:,:,:,:), intent(inout) :: U
    
    U(mu,x(1),x(2),x(3),x(4)) = mat_exp((-dt) * Z(U,x,mu)) *  U(mu,x(1),x(2),x(3),x(4))

  end subroutine gradient_flow

  pure function Z(U,x,mu)
    type(SU2), dimension(:,:,:,:,:), intent(in) :: U
    integer(i4), intent(in) :: x(4), mu
    type(SU2) :: Z
    Z = TA(U(mu,x(1),x(2),x(3),x(4)) * dagger(staples(U,x,mu)))
  end function Z

  pure function TA(W)
    type(SU2), intent(in) :: W
    type(SU2) :: TA, V
    complex(dp) :: trV
    V = W - dagger(W)
    trV = V%matrix(1,1) + V%matrix(2,2)
    TA%matrix = V%matrix/2.0_dp - trV/4 * one%matrix
  end function TA


end module smooth_configurations

  
