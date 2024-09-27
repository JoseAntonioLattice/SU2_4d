module arrays
  use precision
  use datatypes
  implicit none

  type(SU2), allocatable, dimension(:,:,:,:,:) :: U
  real(dp), allocatable, dimension(:) :: beta
  real(dp),allocatable, dimension(:) :: P, Q_den, Poly 
  real(dp) :: avr_P, err_P, avr_Poly, err_Poly
  real(dp) :: avr_Qden, err_Qden,avr_eden, err_eden
  real(dp), allocatable, dimension(:) :: Eden
  real(dp), allocatable, dimension(:,:) :: correlation_polyakov_loop
  real(dp), allocatable, dimension(:) :: avr_cpl, err_cpl
  
contains
  subroutine read_beta()
    character(20) :: file_beta
    integer(i4) :: inunit, stat
    real(dp) :: b

    beta = [real(dp) ::]
    read(*,'(a)') file_beta
    write(*,'(2a)') 'User typed:', trim(file_beta)
    open(newunit = inunit, file = trim(file_beta), iostat = stat, status = 'old')
    
    do
       read(inunit,*,iostat = stat) b 
       if(stat  /= 0) exit
       beta = [beta, b]
    end do
    
    print*, beta
    
  end subroutine read_beta

  subroutine create_beta_array()
    use parameters
    integer(i4) :: i
    real(dp) :: dbeta

    allocate(beta(n_beta))

    dbeta = (b_array(2) - b_array(1))/(n_beta-1)
    beta = [(b_array(1) + (i-1)*dbeta, i = 1, n_beta)]
    
  end subroutine create_beta_array
  
  subroutine reserve_memory
    use parameters
    allocate(U(d,L,L,L,Lt))
    allocate(P(n_measurements))
    allocate(Poly(n_measurements))
    allocate(Q_den(n_measurements))
    allocate(Eden(n_measurements))
    allocate(correlation_polyakov_loop(n_measurements,L/2-1))
    allocate(avr_cpl(L/2-1),err_cpl(L/2-1))
  end subroutine reserve_memory
  
end module arrays
