module arrays
  use precision
  use datatypes
  implicit none

  type(SU2), allocatable, dimension(:,:,:,:,:) :: U
  real(dp), allocatable, dimension(:) :: beta, P
  complex(dp),allocatable, dimension(:) :: Q_den 
  real(dp) :: avr_P, err_P
  complex(dp) :: avr_Qden, err_Qden

contains
  subroutine read_beta()
    character(20) :: file_beta
    integer(i4) ::i, inunit, stat
    real(dp) :: b
    beta = [real(dp) ::]
    read(*,'(a)') file_beta
    write(*,'(2a)') 'User typed:', trim(file_beta)
    
    open(newunit = inunit, file = trim(file_beta), iostat = stat, status = 'old')
    
    print*, stat
    do
       read(inunit,*,iostat = stat) b 
       if(stat  /= 0) exit
       beta = [beta, b]
    end do
    print*, beta
    print*, size(beta)
  end subroutine read_beta
end module arrays
