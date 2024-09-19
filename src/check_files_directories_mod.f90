module check_files_directories_mod

  implicit none
  
contains
  
  subroutine check_directory(directory)
    character(*), intent(in) :: directory
    logical ::  dir_exists

    inquire(file=directory, exist=dir_exists) ! ask wether the directory exists or not
    if(dir_exists .eqv. .false.) call execute_command_line('mkdir '//directory) ! If not it creates it

  end subroutine check_directory

  subroutine check_file(filepath,file_exists)
    character(*), intent(in) :: filepath
    logical, intent(out) :: file_exists

    inquire(file = filepath, exist = file_exists)
    if(file_exists .eqv. .false.)then
       call execute_command_line('touch '//filepath)
       !file_exists = .true.
    end if
    
  end subroutine check_file

  function numbered_file(filename, extension)
    use number2string_mod
    character(*), intent(in) :: filename
    character(*), intent(in) :: extension
    integer :: i
    character(256) :: numbered_file
    logical :: file_exists

    i = 1
    do 
       numbered_file = trim(filename)//"_"//trim(int2str(i))//trim(extension)
       call check_file(trim(numbered_file), file_exists)
       if (file_exists .eqv. .false.) exit
       i = i + 1
    end do
    
  end function numbered_file

  subroutine create_directories(directory_array)
    character(*), dimension(:) :: directory_array
    character(:), allocatable :: dir_name
    integer :: i, n
    
    n = size(directory_array)
    dir_name = ""
    do i = 1, n
       dir_name = trim(dir_name)//trim(directory_array(i))//"/"
       call check_directory(trim(dir_name))
       print*, trim(dir_name)
    end do
    
  end subroutine create_directories


end module check_files_directories_mod
