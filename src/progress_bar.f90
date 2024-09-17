subroutine progress_bar(r)
  use iso_fortran_env, only : outunit => output_unit 
  implicit none
  character :: cr
  real, intent(in) :: r 
  cr = char(13)

  write(outunit, '(a)') 'Progress :'
  write(outunit,'("Progress: ",a,f5.1,"%")',advance = 'no') cr,r*100
  
end subroutine progress_bar
