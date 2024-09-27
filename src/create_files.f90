module create_files
  use check_files_directories_mod
  use number2string_mod
  implicit none

contains

  subroutine create_measurements_file(Lx1,Lt1,beta,algorithm1,equilibrium1)
    use parameters
    integer, intent(in) :: Lx1,Lt1
    real(8), intent(in) :: beta
    character(*), intent(in) :: algorithm1
    logical :: equilibrium1, file_exists, condition
    integer :: i
    character(100), dimension(6) :: directory_array
    
    character(100) :: directory, eq, data_file

    if(equilibrium1 .eqv. .true.) then
       eq = "equilibrium"
    else
       eq = "out_of_equilibrium"
    end if


    directory_array = [ character(100) :: 'data', "Lx="//trim(int2str(Lx1)), &
                                          "Lt="//trim(int2str(Lt1)), trim(eq), &
                                          trim(algorithm1),"beta="//trim(adjustl(real2str(beta))) ]
    
    call create_directories(directory_array)
    directory = ''                                      
    do i = 1, 6
       directory = trim(directory)//trim(directory_array(i))//"/"
    end do
    data_file = trim(directory)//"observables"
    data_file = trim(numbered_file(trim(data_file),".dat"))

    open(unit = 100, file = trim(data_file))
    write(100, nml = input_parameters)
    
  end subroutine create_measurements_file

end module create_files
