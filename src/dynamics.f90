module dynamics

  use precision
  use datatypes
  use local_update_algorithms
  use smooth_configurations
  use observables
  implicit none

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
       q_den(i_sweeps) = topological_charge(U)
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

    write(out_smooth,*) 0, P(0), q_den(0)
    do i_t = 1, N_time
       call sweeps(U,beta,trim(smoothing_method))
       Eden(i_t) = 0.0_dp
       P(i_t) =  plaquette_value(U)
       q_den(i_t) =  topological_charge(U)
       write(out_smooth,*) i_t, P(i_t), q_den(i_t)
       flush(out_smooth)
    end do
    write(out_smooth,'(a)') ' ', ' ', ' '

  end subroutine smooth_configuration
  
  
  subroutine sweeps(U,beta,algorithm)
    use parameters, only : d, L, Lt
    type(SU2), dimension(d,L,L,L,Lt), intent(inout) :: U
    real(dp), intent(in) :: beta
    character(*), intent(in) :: algorithm
    integer(i4) :: x1,x2,x3,x4,mu

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
    case("gradient_flow")
       do x1 = 1, L
          do x2 = 1, L
             do x3 = 1, L
                do x4 = 1, Lt
                   do mu = 1, d
                      call gradient_flow(U,[x1,x2,x3,x4],mu)
                   end do
                end do
             end do
          end do
       end do
       case("cooling")
       do x1 = 1, L
          do x2 = 1, L
             do x3 = 1, L
                do x4 = 1, Lt
                   do mu = 1, d
                      call cooling(U,[x1,x2,x3,x4],mu)
                   end do
                end do
             end do
          end do
       end do
    end select
  end subroutine sweeps
  
end module dynamics
