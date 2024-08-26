program main
  use number2string_mod
  implicit none

  real(8) :: r = 2.0_8, r2 = 20.

  write(*,'(a)') real2str(r)
  write(*,'(a)') real2str(r2)
end program main
