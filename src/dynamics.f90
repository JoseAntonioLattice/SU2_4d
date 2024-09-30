module dynamics

  use precision
  use datatypes
  use local_update_algorithms
  use smooth_configurations
  use observables
  implicit none

  abstract interface
     subroutine lua(u,x,mu,beta)
       use precision
       use datatypes
       type(SU2), dimension(:,:,:,:,:), intent(inout) :: U
       integer(i4), intent(in) :: x(4), mu
       real(dp), intent(in) :: beta
     end subroutine lua
  end interface

  abstract interface
     subroutine lua_nobeta(u,x,mu)
       use precision
       use datatypes
       type(SU2), dimension(:,:,:,:,:), intent(inout) :: U
       integer(i4), intent(in) :: x(4), mu
     end subroutine lua_nobeta
  end interface
  
  abstract interface
     subroutine smooth_configuration_function(U,V,x,mu)
       use precision
       use datatypes
       type(SU2), dimension(:,:,:,:,:), intent(in) :: U
       type(SU2), dimension(:,:,:,:,:), intent(out) :: V
       integer(i4), intent(in) :: x(4), mu
     end subroutine smooth_configuration_function
  end interface
  
contains

  subroutine simulate_equilibrium(U,beta,outunit,P,Poly,correlation_polyakov_loop)
    use starts
    use statistics
    use parameters, only : L,Lt,N_measurements,algorithm
    type(SU2), dimension(:,:,:,:,:), intent(out) :: U
    integer(i4), intent(in) :: outunit
    real(dp), dimension(:), intent(in)  :: beta
    real(dp), dimension(:), intent(out) :: P
    real(dp), dimension(:), intent(out) :: Poly
    real(dp), dimension(:,:), intent(out) :: correlation_polyakov_loop
    integer(i4) :: i, i_b, i_t
    real(dp) :: t1, t2, avr_P, err_P
    
    call cpu_time(t1)
    call progress_bar(0.0)
    call hot_start(U)
    
    i_t = 0
    do i_b = 1, size(beta)
       do i = 1, N_measurements
          i_t = i_t + 1
          call thermalization(U,beta(i_b))
          call measurements(U,P,Poly,correlation_polyakov_loop,beta(i_b),100)
          call progress_bar(i_t/real(N_measurements*size(beta)))
       end do
       call std_err(P,avr_P,err_P)
       write(*,*) char(13),beta(i_b), avr_P, err_P
       write(outunit,*) beta(i_b), avr_P, err_P
       flush(outunit)
    end do
    print*, ''
    call cpu_time(t2)
    print*, "Time: ", t2-t1, "secs"
  end subroutine simulate_equilibrium
  
  subroutine thermalization(U,beta,steps)
    use parameters, only : N_thermalization, algorithm
    type(SU2), dimension(:,:,:,:,:), intent(inout) :: U
    real(dp), intent(in) :: beta
    integer(i4), optional :: steps
    integer(i4) :: i_sweeps, nsteps

    nsteps = N_thermalization
    if( present(steps) ) nsteps = steps
   
    do i_sweeps = 1, nsteps
       call sweeps(U,beta, algorithm)
    end do
    
  end subroutine thermalization

  subroutine measurements(U,P,Poly,corr_poly,beta,unit)
    use parameters, only : L, N_measurements, N_skip,algorithm
    type(SU2), dimension(:,:,:,:,:), intent(inout) :: U
    real(dp), intent(in) :: beta
    real(dp), dimension(N_measurements), intent(out) :: P
    !real(dp), dimension(N_measurements), intent(out) :: Q_den
    !real(dp), dimension(N_measurements), intent(out) :: E_den
    real(dp), dimension(N_measurements,L/2-1), intent(out) :: corr_poly
    !character(*), intent(in) :: definition
    integer(i4), intent(in) :: unit
    real(dp), dimension(N_measurements), intent(out) :: Poly
    integer(i4) :: i_sweeps, i_skip
    real(dp), dimension(L,L,L) :: poly_array
    

    do i_sweeps = 1, N_measurements
       do i_skip = 1, N_skip
          call sweeps(U,beta,algorithm)
       end do
       poly_array = array_polyakov_loop(U)
       Poly(i_sweeps) = sum(poly_array)/L**3
       P(i_sweeps) = plaquette_value(U)
       corr_poly(i_sweeps,:) = correlation(poly_array)
       !write(unit,*) P(i_sweeps), Poly(i_sweeps),corr_poly(i_sweeps,:)
    end do
    
  end subroutine measurements

  subroutine smooth_configuration(U,beta,n_time,smoothing_method,out_smooth)
    type(SU2), dimension(:,:,:,:,:), intent(inout) :: U
    real(dp), intent(in) :: beta
    integer(i4), intent(in) :: n_time
    character(*), intent(in) :: smoothing_method
    integer(i4), intent(in) :: out_smooth 
    integer(i4) :: i_t

    write(out_smooth,*) 0,plaquette_value(U) , topological_charge(U,'plaquette'),&
                        topological_charge(U,'clover'), det(U(1,1,1,1,1)), U(1,1,1,1,1)
    do i_t = 1, N_time
       call sweeps(U,beta,trim(smoothing_method))
       write(out_smooth,*) i_t,plaquette_value(U) , topological_charge(U,'plaquette'),&
                           topological_charge(U,'clover'), det(U(1,1,1,1,1)), U(1,1,1,1,1)
       flush(out_smooth)
    end do
    write(out_smooth,'(a)') ' ', ' ', ' '

  end subroutine smooth_configuration

  subroutine sweeps_lua(U,beta,lua_function)
    use parameters, only : d, L, Lt
    type(SU2), dimension(d,L,L,L,Lt), intent(inout) :: U
    real(dp), intent(in) :: beta
    procedure(lua) :: lua_function
    integer(i4) :: x1,x2,x3,x4,mu
    do x1 = 1, L
       do x2 = 1, L
          do x3 = 1, L
             do x4 = 1, Lt
                do mu = 1, d
                   call lua_function(U,[x1,x2,x3,x4],mu,beta)
                end do
             end do
          end do
       end do
    end do
  end subroutine sweeps_lua

    subroutine sweeps_lua_nobeta(U,lua_function)
    use parameters, only : d, L, Lt
    type(SU2), dimension(d,L,L,L,Lt), intent(inout) :: U
    procedure(lua_nobeta) :: lua_function
    integer(i4) :: x1,x2,x3,x4,mu
    do x1 = 1, L
       do x2 = 1, L
          do x3 = 1, L
             do x4 = 1, Lt
                do mu = 1, d
                   call lua_function(U,[x1,x2,x3,x4],mu)
                end do
             end do
          end do
       end do
    end do
  end subroutine sweeps_lua_nobeta

  subroutine sweeps_scf(U,scf_function)
    use parameters, only : d, L, Lt
    type(SU2), dimension(d,L,L,L,Lt), intent(inout) :: U
    procedure(smooth_configuration_function) :: scf_function
    type(SU2), dimension(d,L,L,L,Lt) :: V
    integer(i4) :: x1,x2,x3,x4,mu
    
    do x1 = 1, L
       do x2 = 1, L
          do x3 = 1, L
             do x4 = 1, Lt
                do mu = 1, d
                   call scf_function(U,V,[x1,x2,x3,x4],mu)
                   U(mu,x1,x2,x3,x4) = V(mu,x1,x2,x3,x4) 
                end do
             end do
          end do
       end do
    end do
    !U = V
  end subroutine sweeps_scf

  subroutine fat_temporal_links(U)
    use parameters, only : d, L, Lt
    type(SU2), dimension(d,L,L,L,Lt), intent(inout) :: U
    type(SU2) :: V
    integer(i4) :: x1,x2,x3,x4,mu

     do x1 = 1, L
       do x2 = 1, L
          do x3 = 1, L
             do x4 = 1, Lt
                V = staples(U,[x1,x2,x3,x4],4)
                U(4,x1,x2,x3,x4) = V/det(V) 
             end do
          end do
       end do
    end do

  end subroutine fat_temporal_links

  subroutine spatial_ape_smearing(U)
    use parameters, only : d, L, Lt
    type(SU2), dimension(d,L,L,L,Lt), intent(inout) :: U
    type(SU2), dimension(d,L,L,L,Lt) :: V
    integer(i4) :: x1, x2, x3, x4, mu
    
    do x1 = 1, L
       do x2 = 1, L
          do x3 = 1, L
             do x4 = 1, Lt
                do mu = 1, d - 1
                   call ape_smearing(U,V,[x1,x2,x3,x4],mu)
                   U(mu,x1,x2,x3,x4) = V(mu,x1,x2,x3,x4)
                end do
             end do
          end do
       end do
    end do
    
  end subroutine spatial_ape_smearing
  
  subroutine sweeps(U,beta,algorithm)
    use parameters, only : d, L, Lt
    type(SU2), dimension(d,L,L,L,Lt), intent(inout) :: U
    real(dp), intent(in) :: beta
    character(*), intent(in) :: algorithm  
   
    select case(algorithm)
    case("heatbath")
       call sweeps_lua(U,beta,heatbath)
    case("metropolis")
       call sweeps_lua(U,beta,metropolis)
    case("overrelaxation")
       call sweeps_lua_nobeta(U,overrelaxation)
    case("gradient_flow")
       call sweeps_scf(U,gradient_flow)
    case("cooling")
       call sweeps_scf(U,cooling)
    case("ape_smearing")
       call sweeps_scf(U,ape_smearing)
    end select
    
  end subroutine sweeps
  
end module dynamics
