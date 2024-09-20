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
  
  subroutine thermalization(U,beta)
    use parameters, only : N_thermalization, algorithm
    type(SU2), dimension(:,:,:,:,:), intent(inout) :: U
    real(dp), intent(in) :: beta
    integer(i4) :: i_sweeps

    do i_sweeps = 1, N_thermalization
       call sweeps(U,beta, algorithm)
    end do
    
  end subroutine thermalization

  subroutine measurements(U,beta,P,Q_den, E_den,definition)
    use parameters, only : N_measurements, N_skip,algorithm
    type(SU2), dimension(:,:,:,:,:), intent(inout) :: U
    real(dp), intent(in) :: beta
    real(dp), dimension(N_measurements), intent(out) :: P
    real(dp), dimension(N_measurements), intent(out) :: Q_den
    real(dp), dimension(N_measurements), intent(out) :: E_den
    character(*), intent(in) :: definition
    integer(i4) :: i_sweeps, i_skip

    do i_sweeps = 1, N_measurements
       do i_skip = 1, N_skip
          call sweeps(U,beta,algorithm)
       end do
       P(i_sweeps) = plaquette_value(U)
       q_den(i_sweeps) = topological_charge(U,'plaquette')
       E_den(i_sweeps) = E(U,definition)
       write(100,*) P(i_sweeps), q_den(i_sweeps), E_den(i_sweeps)
    end do
    
  end subroutine measurements

  subroutine smooth_configuration(U,beta,n_time,eden,P,q_den,smoothing_method,out_smooth)
    type(SU2), dimension(:,:,:,:,:), intent(inout) :: U
    real(dp), intent(in) :: beta
    integer(i4), intent(in) :: n_time
    real(dp), dimension(0:), intent(inout) :: P
    real(dp), dimension(0:), intent(inout) :: Q_den
    real(dp), dimension(0:), intent(inout) :: Eden
    character(*), intent(in) :: smoothing_method
    integer(i4), intent(in) :: out_smooth 
    integer(i4) :: i_t

    write(out_smooth,*) 0, P(0), q_den(0), topological_charge(U,'clover'), det(U(1,1,1,1,1)), U(1,1,1,1,1)
    do i_t = 1, N_time
       call sweeps(U,beta,trim(smoothing_method))
       !print*, det(U(1,1,1,1,1))
       Eden(i_t) = 0.0_dp
       P(i_t) =  plaquette_value(U)
       q_den(i_t) =  topological_charge(U,'plaquette')
       write(out_smooth,*) i_t, P(i_t), q_den(i_t), topological_charge(U,'clover'), det(U(1,1,1,1,1)), U(1,1,1,1,1)
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
