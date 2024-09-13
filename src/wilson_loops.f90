module wilson_loops
  use precision
  use pbc, only : ip, im
  use datatypes
  implicit none

contains

  pure function plaquette(U,x,mu,nu)
    type(SU2), dimension(:,:,:,:,:), intent(in) :: U
    integer(i4), intent(in) :: x(4), mu, nu
    integer(i4) :: x2(4),x3(4)
    type(SU2) :: plaquette

    x2 = ip(x,mu)
    x3 = ip(x,nu)

    plaquette = U(mu, x(1), x(2), x(3), x(4)) * &
                U(nu,x2(1),x2(2),x2(3),x2(4)) * &
         dagger(U(mu,x3(1),x3(2),x3(3),x3(4)))* &
         dagger(U(nu, x(1), x(2), x(3), x(4))) 
  end function plaquette

  pure function clover(U,x,mu,nu)
    type(SU2) :: clover
    type(SU2), dimension(:,:,:,:,:), intent(in) :: U
    integer(i4), intent(in) :: x(4), mu, nu
    integer(i4), dimension(4) :: x2, x3,x4,x5,x6,x7,x8 

    x2 = ip(x,mu)  ! x + mu
    x3 = ip(x,nu)  ! x + nu
    x4 = im(x,mu)  ! x - mu
    x5 = im(x,nu)  ! x - nu
    x6 = ip(x4,nu) ! x - mu + nu
    x7 = im(x4,nu) ! x - mu - nu
    x8 = im(x2,nu) ! x + mu - nu 
    clover = plaquette(U,x,mu,nu) + &
         U(nu,x(1),x(2),x(3),x(4)) * dagger(U(mu,x6(1),x6(2),x6(3),x6(4))) * &
         dagger(U(nu,x4(1),x4(2),x4(3),x4(4))) * U(mu,x4(1),x4(2),x4(3),x4(4)) + &
         dagger(U(mu,x4(1),x4(2),x4(3),x4(4))) * dagger(U(nu,x7(1),x7(2),x7(3),x7(4))) * &
         U(mu,x7(1),x7(2),x7(3),x7(4)) *  U(nu,x5(1),x5(2),x5(3),x5(4)) + &
         dagger(U(nu,x5(1),x5(2),x5(3),x5(4))) * U(mu,x5(1),x5(2),x5(3),x5(4)) * &
         U(nu,x8(1),x8(2),x8(3),x8(4)) * dagger(U(mu,x(1),x(2),x(3),x(4)))

  end function clover

  pure function staples(U,x,mu)
    use parameters, only: d
    type(SU2), dimension(:,:,:,:,:), intent(in) :: U
    integer(i4), intent(in) :: x(4), mu
    integer(i4) :: nu
    integer(i4), dimension(4) :: x2,x3,x4,x5
    type(SU2) :: staples

    staples%matrix = 0.0_dp
    x3 = ip(x,mu) ! x + mu
    do nu = 1, d
       if( nu == mu) cycle
       x2 = ip( x,nu) ! x + nu
       x4 = im( x,nu) ! x - nu
       x5 = im(x3,nu) ! x + mu + nu
       staples = staples + U(nu, x(1), x(2), x(3), x(4)) *&
                           U(mu,x2(1),x2(2),x2(3),x2(4)) *&
                    dagger(U(nu,x3(1),x3(2),x3(3),x3(4)))+&
                    dagger(U(nu,x4(1),x4(2),x4(3),x4(4)))*&
                           U(mu,x4(1),x4(2),x4(3),x4(4)) *&
                           U(nu,x5(1),x5(2),x5(3),x5(4))
    end do
  end function staples
    
end module wilson_loops
