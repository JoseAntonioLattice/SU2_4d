module dynamics

  use data_types
  use pbc 
  use iso_fortran_env, only : dp => real64, i4 => int32
  implicit none

  real(dp), parameter :: pi = acos(-1.0_dp)
  
contains

  subroutine set_memory(U,beta,L,Lt,N_beta,bi,bf,plqaction,n_measurements)
    type(SU2), allocatable, dimension(:,:,:,:,:) :: U
    real(dp), allocatable, dimension(:) :: beta
    integer(i4), intent(in) :: L,Lt,N_beta,n_measurements
    real(dp), intent(in) :: bi, bf
    real(dp), allocatable, dimension(:) :: plqaction
   
    
    integer(i4) :: i

    call set_periodic_bounds(L,Lt)
    allocate(U(4,L,L,L,Lt))
    allocate(beta(N_beta))
    allocate(plqaction(n_measurements))
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

  function staples(U,x,mu) result(A)

    type(SU2), dimension(:,:,:,:,:), intent(in) :: U
    integer(i4), intent(in) :: x(4), mu
    integer(i4) :: nu
    integer(i4), parameter :: d = 4
    type(SU2) :: A

    integer(i4), dimension(4) :: x2,x3,x4,x5

    A%matrix = (0.0_dp,0.0_dp)
    x3 = ip_func(x,mu) ! x + mu
    do nu = 1, d
       if(nu .ne. mu)then

          x2 = ip_func(x,nu) ! x + nu
          x4 = im_func(x,nu) ! x - nu
          x5 = ip_func(x4,mu)! x - nu + mu

          A = A + U(nu,x(1),x(2),x(3),x(4))*U(mu,x2(1),x2(2),x2(3),x2(4))!*dagger(U(nu,x3(1),x3(2),x3(3),x3(4))!+ &
             !dagger(U(nu,x4(1),x4(2),x4(3),x4(4))*U(mu,x4(1),x4(2),x4(3),x4(4))*U(nu,x5(1),x5(2),x5(3),x5(4))
       end if
    end do
  end function staples


  subroutine heatbath(U,x,mu,beta)

    type(SU2), dimension(:,:,:,:,:), intent(inout) :: U
    integer(i4), intent(in) :: mu, x(4)
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
    n(2) = sqrt(1.0_dp - n(1)**2) * cos(2*pi*W(2))
    n(3) = sqrt(1.0_dp - n(1)**2) * sin(2*pi*W(2))
    
    n(1:3) = sqrt(1.0_dp - n(0)**2) * n(1:3)
    
    a = cmplx(n(0),n(1),dp)
    b = cmplx(n(2),n(3),dp)
    
    Xmat%matrix = reshape([a,-conjg(b),b,conjg(a)],[2,2])
    
    U(mu,x(1),x(2),x(3),x(4)) = Xmat * V
    
  end subroutine heatbath

  subroutine sweeps(u,beta)

    type(SU2), dimension(:,:,:,:,:), intent(inout) :: u
    real(dp), intent(in) :: beta

    integer(i4) :: x,y,z,t,mu
    integer(i4) :: L,Lt

    L  = size(u(1,:,1,1,1))
    Lt = size(u(1,1,1,1,:))
    do x = 1, L
       do y = 1, L
          do z = 1, L
             do t = 1, Lt
                do mu = 1, 4
                   call heatbath(U,[x,y,z,t],mu,beta)
                end do
             end do
          end do
       end do
    end do

  end subroutine sweeps

  subroutine initialization(u,plqaction,beta,N_thermalization,N_measurements, N_skip)

    type(SU2), dimension(:,:,:,:,:), intent(inout) :: u
    real(dp), dimension(:), intent(out) :: plqaction
    integer(i4), intent(in) :: N_thermalization, N_measurements, N_skip
    real(dp), intent(in) :: beta

    integer(i4) :: i_skip, i_sweeps

    call thermalization(u, N_thermalization, beta)
    
    do i_sweeps = 1, N_measurements
       do i_skip = 1, n_skip
          call sweeps(u,beta)
       end do
       plqaction(i_sweeps) = action(u)
    end do

  end subroutine initialization

  subroutine thermalization(u,N_thermalization,beta)

    type(SU2), dimension(:,:,:,:,:), intent(inout) :: u
    integer(i4) :: N_thermalization
    real(dp) :: beta

    integer(i4) :: i_sweeps
    
    do i_sweeps = 1, N_thermalization
       call sweeps(u,beta)
    end do
    
  end subroutine thermalization

function action(u)

    real(dp) :: action

    type(SU2), dimension(:,:,:,:,:), intent(in) :: u

    integer(i4) :: x,y,z,t,mu,nu
    integer(i4) :: L,Lt

    L = size(U(1,:,1,1,1))
    Lt = size(U(1,1,1,1,:))
    action = 0.0_dp

    do x = 1, L
       do y = 1, L
          do z = 1, L
             do t = 1, Lt
                do mu = 1, 3
                   do nu = mu+1,4
                      action = action + real(tr(plaquette(u,[x,y,z,t],mu,nu)))
                   end do
                end do
             end do
          end do
       end do
    end do
    action = action / (L**3*Lt)

  end function action

  function plaquette(u,x,mu,nu)

    type(SU2) :: plaquette

    type(SU2), dimension(:,:,:,:,:), intent(in) :: u
    integer(i4), intent(in) :: x(4),mu,nu
    integer(i4), dimension(4) :: x2, x3


    x2 = ip_func(x,mu)
    x3 = ip_func(x,nu)

    plaquette = U(mu,x(1),x(2),x(3),x(4)) * U(nu,x2(1),x2(2),x2(3),x(4)) * &
          dagger(U(mu,x3(1),x3(2),x3(3),x(4))) * dagger(U(nu,x(1),x(2),x(3),x(4)))

  end function plaquette

end module dynamics
