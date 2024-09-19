module configurations
  use datatypes
  use parameters, only : d, L, Lt
  use number2string_mod
  use check_files_directories_mod
  implicit none

contains

  subroutine save_configuration(U,beta)
    type(SU2), dimension(d,L,L,L,Lt), intent(in)  :: U
    real(dp), intent(in) :: beta
    character(256) :: configurations_file
    integer(i4) :: i, configurations_unit
    character(100), dimension(4) :: directory_array

    directory_array = [character(100):: "configurations", "Lx="//trim(int2str(L)), "Lt="//trim(int2str(Lt)), &
                       "beta="//trim(real2str(beta))]
    call create_directories(directory_array)
    configurations_file = ''
    do i = 1, 4
       configurations_file = trim(configurations_file)//trim(directory_array(i))//"/"
    end do
    configurations_file = trim(configurations_file)//"configuration"
    configurations_file = trim(numbered_file(trim(configurations_file),".conf")) 
    print '(a,i0,a)', "The size of the file is ",4*d*L**3*Lt*16, " bytes."
    open(newunit = configurations_unit, file = trim(configurations_file), form = "unformatted", access = "sequential")
    write(configurations_unit) U
    close(configurations_unit)
  end subroutine save_configuration

  subroutine read_configuration(U,beta,i)
    type(SU2), dimension(d,L,L,L,Lt), intent(out) :: U
    real(dp), intent(in) :: beta
    integer(i4), intent(in) :: i
    character(256) :: configurations_file
    integer(i4) :: j, configurations_unit
    character(100), dimension(4) :: directory_array
    
    directory_array = [character(100):: "configurations", "Lx="//trim(int2str(L)), "Lt="//trim(int2str(Lt)), &
                       "beta="//trim(real2str(beta))]
    configurations_file = ''
    do j = 1, 4
       configurations_file = trim(configurations_file)//trim(directory_array(j))//"/"
    end do
    configurations_file = trim(configurations_file)//"configuration_"//trim(int2str(i))//".conf"
    print*, trim(configurations_file)
    open(newunit = configurations_unit, file = trim(configurations_file), form = "unformatted", access = "sequential")
    read(configurations_unit) U
    close(configurations_unit)
  end subroutine read_configuration

end module configurations
