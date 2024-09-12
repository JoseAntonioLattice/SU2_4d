module starts

  use datatypes
  implicit none

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
end module starts
