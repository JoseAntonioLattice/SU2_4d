module number2string_mod

    use iso_fortran_env, only : dp => real64, i4 => int32

    implicit none

    private

    public :: int2str, real2str

contains

    character(len=20) function int2str(k)
      !"Convert an integer to string."
      integer(i4), intent(in) :: k
      write (int2str, *) k
      int2str = adjustl(int2str)
    end function int2str

    character(20) function real2str(r)
      real(dp), intent(in) :: r
      
      write(real2str,"(f10.4)") r
      real2str = adjustl(real2str)
      
    end function real2str

end module number2string_mod
