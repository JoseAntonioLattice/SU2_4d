module dynamics

  use precision
  use datatypes
  use pbc, only : ip, im
  implicit none

  type(SU2) :: one
  real(dp), parameter :: pi = acos(-1.0_dp)
contains

  ! STARTS
  subroutine hot_start(U)
    use parameters, only : d, L, Lt
    type(SU2), intent(out) :: U(d,L,L,L,Lt)
    integer(i4) :: x1,x2,x3,x4,mu

    do x1 = 1, L
       do x2 = 1, L
          do x3 = 1, L
             do x4 = 1, Lt
                do mu = 1, d
                   U(mu,x1,x2,x3,x4) = SU2_ran()
                end do
             end do
          end do
       end do
    end do
    
  end subroutine hot_start

  subroutine cold_start(U)
    use parameters, only : d, L, Lt
    type(SU2), intent(out) :: U(d,L,L,L,Lt)
    U = one
  end subroutine cold_start

  subroutine create_one()
    one = SU2_matrix((1.0_dp,0.0_dp),(0.0_dp,0.0_dp))
  end subroutine create_one
  
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

  function SU2_matrix(a,b)
    type(SU2) :: SU2_matrix
    complex(dp) :: a,b

    SU2_matrix%matrix(1,1) = a
    SU2_matrix%matrix(1,2) = b
    SU2_matrix%matrix(2,1) = -conjg(b)
    SU2_matrix%matrix(2,2) =  conjg(a)
  end function SU2_matrix

  function plaquette(U,x,mu,nu)
    type(SU2), dimension(:,:,:,:,:), intent(in) :: U
    integer(i4), intent(in) :: x(4), mu
    integer(i4) :: nu,x2(4),x3(4)
    type(SU2) :: plaquette

    x2 = ip(x,mu)
    x3 = ip(x,nu)

    plaquette = U(mu, x(1), x(2), x(3), x(4)) * &
                U(nu,x2(1),x2(2),x2(3),x2(4)) * &
         dagger(U(mu,x3(1),x3(2),x3(3),x3(4)))* &
         dagger(U(nu, x(1), x(2), x(3), x(4))) 
  end function plaquette
  
  function plaquette_value(U) result(P)
    use parameters, only : d, N, L, Lt
    type(SU2), dimension(d,L,L,L,Lt), intent(in) :: U
    real(dp) :: P
    integer(i4) :: vol
    integer(i4), parameter :: planes = d*(d-1)/2 
    integer(i4) :: x1,x2,x3,x4,mu,nu

    P = 0.0_dp
    do x1 = 1, L
       do x2 = 1, L
          do x3 =1, L
             do x4 = 1, Lt
                do mu = 1, d - 1
                   do nu = mu + 1, d
                      P = P + real(tr(plaquette(U,[x1,x2,x3,x4],mu,nu)))
                   end do
                end do
             end do
          end do
       end do
    end do
    vol = L**3*Lt
    P = P/(N*planes*vol)
  end function plaquette_value
  
  function staples(U,x,mu)
    use parameters, only: d
    type(SU2), dimension(:,:,:,:,:), intent(in) :: U
    integer(i4), intent(in) :: x(4), mu
    integer(i4) :: nu
    integer(i4), dimension(4) :: x2,x3,x4,x5
    type(SU2) :: staples

    staples%matrix = 0.0_dp
    x3 = ip(x,mu) ! x + mu
    do nu = 1, d
       if( nu .ne. mu) then
          x2 = ip( x,nu) ! x + nu
          x4 = im( x,nu) ! x - nu
          x5 = im(x3,nu) ! x + mu + nu
          staples = staples + U(nu, x(1), x(2), x(3), x(4)) *&
                              U(mu,x2(1),x2(2),x2(3),x2(4)) *&
                       dagger(U(nu,x3(1),x3(2),x3(3),x3(4)))+&
                       dagger(U(nu,x4(1),x4(2),x4(3),x4(4)))*&
                              U(mu,x4(1),x4(2),x4(3),x4(4)) *&
                              U(nu,x5(1),x5(2),x5(3),x5(4))
       end if
    end do
  end function staples

  function DS(U,Up,x,mu,beta)
    use parameters, only : N
    type(SU2), dimension(:,:,:,:,:), intent(in) :: U
    type(SU2), intent(in) :: Up
    integer(i4), intent(in) :: x(4), mu
    real(dp), intent(in) :: beta
    real(dp) :: DS

    DS = -(beta/N)*real( tr( (Up - U(mu,x(1),x(2),x(3),x(4))) * dagger(staples(U,[x(1),x(2),x(3),x(4)],mu)) ))

  end function DS

  subroutine thermalization(U,beta)
    use parameters, only : N_thermalization
    type(SU2), dimension(:,:,:,:,:), intent(inout) :: U
    real(dp), intent(in) :: beta
    integer(i4) :: i_sweeps

    do i_sweeps = 1, N_thermalization
       call sweeps(U,beta)
    end do
  end subroutine thermalization

  subroutine measurements(U,beta,P)
    use parameters, only : N_measurements, N_skip
    type(SU2), dimension(:,:,:,:,:), intent(inout) :: U
    real(dp), intent(in) :: beta
    real(dp), dimension(N_measurements), intent(out) :: P
    integer(i4) :: i_sweeps, i_skip

    do i_sweeps = 1, N_measurements
       do i_skip = 1, N_skip
          call sweeps(U,beta)
       end do
       P(i_sweeps) = plaquette_value(U)
       write(100,*) P(i_sweeps)
    end do
  end subroutine measurements
  
  subroutine sweeps(U,beta)
    use parameters, only : d, L, Lt, algorithm
    type(SU2), dimension(d,L,L,L,Lt), intent(inout) :: U
    real(dp), intent(in) :: beta
    integer(i4) :: x1,x2,x3,x4,mu
    type(SU2) :: Up

    select case(algorithm)
    case("heatbath")
       do x1 = 1, L
          do x2 = 1, L
             do x3 = 1, L
                do x4 = 1, Lt
                   do mu = 1, d
                      call heatbath(U,[x1,x2,x3,x4],mu,beta)
                   end do
                end do
             end do
          end do
       end do
    case("metropolis")
       do x1 = 1, L
          do x2 = 1, L
             do x3 = 1, L
                do x4 = 1, Lt
                   do mu = 1, d
                      call metropolis(U,[x1,x2,x3,x4],mu,beta)
                   end do
                end do
             end do
          end do
       end do
       case("overrelaxation")
       do x1 = 1, L
          do x2 = 1, L
             do x3 = 1, L
                do x4 = 1, Lt
                   do mu = 1, d
                      call overrelaxation(U,[x1,x2,x3,x4],mu)
                   end do
                end do
             end do
          end do
       end do 
    end select
  end subroutine sweeps
  
  subroutine metropolis(U,x,mu,beta)
    type(SU2), dimension(:,:,:,:,:), intent(inout) :: U
    integer(i4) :: x(4), mu
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
    use parameters, only : L, Lt
    type(SU2), dimension(4,L,L,L,Lt), intent(inout) :: U
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
    integer(i4), intent(in) :: mu, x(4)
    type(SU2) :: sigma,V
    real(dp) :: alpha

    sigma = staples(U,x,mu)
    alpha = sqrt(det(sigma))
    V = sigma/alpha
    U(mu,x(1),x(2),x(3),x(4)) = V * dagger(U(mu,x(1),x(2),x(3),x(4))) * V
    !U(mu,x(1),x(2),x(3),x(4)) = tr(U(mu,x(1),x(2),x(3),x(4))*dagger(V))* V - U(mu,x(1),x(2),x(3),x(4))
  end subroutine overrelaxation

  function TA(W)
    type(SU2), intent(in) :: W
    type(SU2) :: TA
    TA = (W - dagger(W))/2.0_dp - tr(W - dagger(W))/(4.0_dp) * one
  end function TA

  function Z(U,x,mu)
    integer(i4), intent(in) :: x(4), mu
    type(SU2) :: Z
    Z = -1.0_dp * TA(U(mu,x(1),x(2),x(3),x(4)) * dagger(staples(U,[x(1),x(2),x(3),x(4)],mu)))
  end function Z

  subroutine Wilson_flow(U,x,mu)
    integer(i4), intent(inout) :: x(4), mu
    type(SU2), intent(in) :: U
    type(SU2) :: V
    integer(i4), parameter :: dt = 0.1
    integer(i4) :: i
    
    V = U(mu,x(1),x(2),x(3),x(4))
    do i = 1, 100
       V = mat_exp(dt * Z(U,x,mu)) * V
    end do
    U(mu,x(1),x(2),x(3),x(4)) = V
  end subroutine Wilson_flow

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
    B%matrix = one
    B%matrix(1,1) = exp(eigenv(1))
    B%matrix(2,2) = exp(eigenv(2))
    B%matrix(3,3) = exp(eigenv(3))

    res = C*B*inv(C)

  end function mat_exp

  function Q(U,x,mu,nu)
    type(SU2) :: Q
    type(SU2), dimension(:,:,:,:,:), intent(in) :: U
    integer(i4), intent(in) :: x(4), mu
    integer(i4), dimension(4) :: x2, x3,x4,x5,x6,x7,x8 

    x2 = ip(x,mu)  ! x + mu
    x3 = ip(x,nu)  ! x + nu
    x4 = im(x,mu)  ! x - mu
    x5 = im(x,nu)  ! x - nu
    x6 = ip(x4,nu) ! x - mu + nu
    x7 = im(x4,nu) ! x - mu - nu
    x8 = im(x2,nu) ! x + mu - nu 
    Q = plaquette(U,x,mu,nu) + &
        U(nu,x(1),x(2),x(3),x(4)) * dagger(U(mu,x6(1),x6(2),x6(3),x6(4))) * dagger(U(nu,x4(1),x4(2),x4(3),x4(4))) * U(mu,x4(1),x4(2),x4(3),x4(4)) + &
        dagger(U(mu,x4(1),x4(2),x4(3),x4(4))) * dagger(U(nu,x7(1),x7(2),x7(3),x7(4))) *  U(mu,x7(1),x7(2),x7(3),x7(4)) *  U(nu,x5(1),x5(2),x5(3),x5(4)) + &
        dagger(U(nu,x5(1),x5(2),x5(3),x5(4))) * U(mu,x5(1),x5(2),x5(3),x5(4)) *  U(nu,x8(1),x8(2),x8(3),x8(4)) * dagger(U(mu,x(1),x(2),x(3),x(4)))
         
  end function Q

  
  function F(U,x,mu,nu)
    type(SU2) :: F, top
    type(SU2), dimension(:,:,:,:,:), intent(in) :: U
    integer(i4), intent(in) :: x(4), mu

    top = Q(U,x,mu,nu) 
    F = (top - dagger(top))/8.0_dp
         
  end function F

  function dual(A,mu,nu)
    type(SU2), intent(in) :: A
    type(SU2) :: dual
    integer(i4), intent(in) :: mu, nu

    dual%matrix = (0.0_dp,0.0_dp)
    
  end function dual
  
end module dynamics
