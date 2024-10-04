module observables
  use precision
  use datatypes, id => levi_civita_indices 
  use wilson_loops

  implicit none

  abstract interface
     pure function top_char_den_function(U,x)
       use precision
       use datatypes
       type(SU2), dimension(:,:,:,:,:), intent(in) :: U
       integer(i4), dimension(4), intent(in) :: x
       real(dp) :: top_char_den_function
     end function top_char_den_function
  end interface
  
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
          do x3 = 1, L
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
          do x3 = 1, L
             do x4 = 1, Lt
                do mu = 1, d - 1
                   do nu = mu + 1, d
                      energy_density = energy_density + tr(one - plaquette(U,[x1,x2,x3,x4],mu,nu))
                   end do
                end do
             end do
          end do
       end do
    end do
    vol = L**3*Lt
    energy_density = energy_density/(3*vol)
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
  
  pure function E(U)
    use parameters, only: d,L,Lt
    real(dp) :: E
    type(SU2), intent(in), dimension(d,L,L,L,Lt) :: U
    !character(*), intent(in) :: definition
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
                      E = E + tr(F_clover(U,x,mu,nu)*F_clover(U,x,nu,mu))
                   end do
                end do
             end do
          end do
       end do
    end do
    Vol = L**3*Lt
    E = -E / 64
    E = E/(2*Vol)
  end function E

  pure function F_clover(U,x,mu,nu)
    type(SU2), dimension(:,:,:,:,:), intent(in) :: U
    integer(i4), dimension(4), intent(in) :: x
    integer(i4), intent(in) :: mu, nu
    type(SU2) :: F_clover
    
    F_clover = clover(U,x,mu,nu) -  clover(U,x,nu,mu)

  end function F_clover

  
  pure function top_den_clover(U,x)
    type(SU2), dimension(:,:,:,:,:), intent(in) :: U
    integer(i4), dimension(4), intent(in) :: x
    real(dp) :: top_den_clover
    type(SU2) :: clov1, clov2, clov3
    integer(i4) :: mu, nu, rho, sigma

    mu = 1; nu = 2; rho = 3; sigma = 4
    
    clov1 = clover(U,x,rho,sigma)
    clov2 = clover(U,x,nu ,sigma)
    clov3 = clover(U,x,nu ,rho  )
    top_den_clover = tr( clover(U,x,mu,nu   ) * (clov1 - dagger(clov1)) ) &
                    -tr( clover(U,x,mu,rho  ) * (clov2 - dagger(clov2)) ) &
                    +tr( clover(U,x,mu,sigma) * (clov3 - dagger(clov3)) )  
            
  end function top_den_clover

  pure function top_den_plaquette(U,x)
    type(SU2), dimension(:,:,:,:,:), intent(in) :: U
    integer(i4), dimension(4), intent(in) :: x
    real(dp) :: top_den_plaquette
    type(SU2) :: clov1, clov2, clov3
    integer(i4) :: i, mu, nu, rho, sigma

    top_den_plaquette = 0.0_dp
    do i = 1, 24
       mu    = id(i,1)
       nu    = id(i,2)
       rho   = id(i,3)
       sigma = id(i,4)
       top_den_plaquette = top_den_plaquette + levi_civita(mu,nu,rho,sigma) * &
                           tr(plaquette(U,x,mu,nu)*plaquette(U,x,rho,sigma))
    end do
  end function top_den_plaquette
  
  pure function topological_charge(U,definition)
    use parameters, only: L,Lt,d
    type(SU2), dimension(d,L,L,L,Lt), intent(in) :: U
    character(*), intent(in) :: definition
    integer(i4) :: x, y, z, t
    real(dp) :: topological_charge
    real(dp), parameter :: pi = acos(-1.0_dp)

    select case(definition)
    case('plaquette')
       topological_charge = -sum_topological_charge_density(U,top_den_plaquette)/(32*pi**2)
    case('clover')
       topological_charge = -sum_topological_charge_density(U,top_den_clover)/(128*pi**2)
    end select
  end function topological_charge

  pure function sum_topological_charge_density(U,tcdf) result(top)
    use parameters, only: L,Lt,d
    type(SU2), dimension(d,L,L,L,Lt), intent(in) :: U
    procedure(top_char_den_function) :: tcdf
    real(dp) :: top
    integer(i4) :: x, y, z, t
    real(dp), dimension(L,L,L,Lt) :: top_char

    do x = 1, L
       do y = 1, L
          do z = 1, L
             do t = 1, Lt
                top_char(x,y,z,t) = tcdf(U,[x,y,z,t])
             end do
          end do
       end do
    end do
    top = sum(top_char)
    
  end function sum_topological_charge_density

  pure function array_polyakov_loop(U)
    use parameters, only : d, L, Lt
    real(dp), dimension(L,L,L) :: array_polyakov_loop 
    type(SU2), dimension(d,L,L,L,Lt), intent(in) :: U
    integer(i4) :: x, y, z

    do x = 1, L
       do y = 1, L
          do z = 1, L
             array_polyakov_loop(x,y,z) = polyakov_loop(U,[x,y,z])
          end do
       end do
    end do
    
  end function array_polyakov_loop

  pure function correlation(array)
    use parameters, only : L
    real(dp), dimension(L/2-1) :: correlation
    real(dp), dimension(L,L,L), intent(in) :: array 
    real(dp), dimension(L,L,L) :: array_corr
    integer(i4) :: x, y, z, t, xp, yp, zp 

    do t = 1, L/2 - 1
       do x = 1, L
          xp = mod(x+t,L); if(xp == 0) xp = L
          do y = 1, L
             yp = mod(y+t,L); if(yp == 0) yp = L
             do z = 1, L
                zp = mod(z+t,L); if(zp == 0) zp = L
                array_corr(x,y,z) = array(x,y,z) * (array(xp,y,z) + array(x,yp,z) + array(x,y,zp))
             end do
          end do
       end do
       correlation(t) = sum(array_corr)/(3*L**3)
    end do
    
  end function correlation
  
end module observables
