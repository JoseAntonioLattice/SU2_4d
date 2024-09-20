program create_dir

  implicit none
  character(10), dimension(4) :: dir_arr

  dir_arr = [character(10) :: "hola","beibi","te","amo"]
  call create_directories(dir_arr)
  
contains

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
end program create_dir
