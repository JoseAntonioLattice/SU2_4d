module observables
  use precision
  use datatypes, id => levi_civita_indices 
  use wilson_loops
  implicit none
contains

  pure function plaquette_value(U) result(P)
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

  pure function energy_density(U) 
    use parameters, only : d, L, Lt
    type(SU2), dimension(d,L,L,L,Lt), intent(in) :: U
    real(dp) :: energy_density
    integer(i4) :: vol
    integer(i4) :: x1,x2,x3,x4,mu,nu

    energy_density = 0.0_dp
    
    do x1 = 1, L
       do x2 = 1, L
          do x3 =1, L
             do x4 = 1, Lt
                do mu = 1, d - 1
                   do nu = mu + 1, d
                      energy_density = energy_density + real(tr(one - plaquette(U,[x1,x2,x3,x4],mu,nu)))
                   end do
                end do
             end do
          end do
       end do
    end do
    vol = L**3*Lt
    energy_density = 2*energy_density/vol
  end function energy_density
  
  pure function DS(U,Up,x,mu,beta)
    use parameters, only : N
    type(SU2), dimension(:,:,:,:,:), intent(in) :: U
    type(SU2), intent(in) :: Up
    integer(i4), intent(in) :: x(4), mu
    real(dp), intent(in) :: beta
    real(dp) :: DS

    DS = -(beta/N)*real( tr( (Up - U(mu,x(1),x(2),x(3),x(4))) * dagger(staples(U,[x(1),x(2),x(3),x(4)],mu)) ))

  end function DS

  pure function F(U,x,mu,nu,definition)
    !real(dp), dimension(2,2) :: F
    type(SU2) :: F, top
    type(SU2), dimension(:,:,:,:,:), intent(in) :: U
    integer(i4), intent(in) :: x(4), mu, nu
    character(*), intent(in) :: definition
    !complex(dp) :: i = (0.0_dp,1.0_dp)
    
    select case(definition)
    case('plaquette')
       F = plaquette(U,x,mu,nu)
       !F = aimag(top%matrix)
       
    case('clover')
       top = clover(U,x,mu,nu)
       !F = aimag(top%matrix)/4.0
       F = (top - dagger(top))/8.0_dp
    end select
  end function F
  
  pure function E(U,definition)
    use parameters, only: d,L,Lt
    real(dp) :: E
    type(SU2), intent(in), dimension(d,L,L,L,Lt) :: U
    character(*), intent(in) :: definition
    integer(i4) :: x1,x2,x3,x4,x(4), mu, nu
    integer(i4) :: Vol

    E = 0.0_dp
    do x1 = 1, L
       do x2 = 1, L
          do x3 = 1, L
             do x4 = 1, Lt
                x = [x1,x2,x3,x4]
                do mu = 1, d
                   do nu = 1, d
                      E = E + tr(F(U,x,mu,nu,definition)*F(U,x,mu,nu,definition))
                   end do
                end do
             end do
          end do
       end do
    end do
    Vol = L**3*Lt
    E = -E/(2*Vol)
  end function E
  
  pure function top_den(U,x)
    type(SU2), dimension(:,:,:,:,:), intent(in) :: U
    real(dp) :: top_den
    integer(i4), dimension(4), intent(in) :: x
    type(SU2) :: clov1, clov2, clov3
    integer(i4), parameter :: mu = 1, nu = 2, rho = 3, sigma = 4
    
    clov1 = clover(U,x,rho,sigma)
    clov2 = clover(U,x,nu ,sigma)
    clov3 = clover(U,x,nu ,rho  )
    top_den = real(tr( clover(U,x,mu,nu )   * (clov1 - dagger(clov1)) )) &
             -real(tr( clover(U,x,mu,rho)   * (clov2 - dagger(clov2)) )) &
             +real(tr( clover(U,x,mu,sigma) * (clov3 - dagger(clov3)) ))  

  end function top_den
  
  pure function topological_charge(U)
    use parameters, only: L,Lt,d
    type(SU2), dimension(d,L,L,L,Lt), intent(in) :: U
    integer(i4) :: x, y, z, t
    real(dp) :: topological_charge
    real(dp), parameter :: pi = acos(-1.0_dp)
    topological_charge = 0.0_dp
    do x = 1, L
       do y = 1, L
          do z = 1, L
             do t = 1, Lt
                topological_charge = topological_charge + top_den(U,[x,y,z,t])
             end do
          end do
       end do
    end do
    topological_charge = -topological_charge/(128*pi**2)                  
  end function topological_charge

end module observables
