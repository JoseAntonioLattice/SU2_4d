subroutine progress_bar(r)
  use iso_fortran_env, only : outunit => output_unit 
  implicit none
  character :: cr
  real, intent(in) :: r 
  cr = char(13)

  write(outunit,'("Progress: ",2a,f5.1,"%")',advance = 'no') cr,"Progress: ",r*100
  
end subroutine progress_bar
