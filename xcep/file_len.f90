  ! 
  ! check file length  
  !
  integer function file_len(filename)

    implicit none 

    character(len=*) :: filename
    character(len=100) :: ss
   ! integer :: file_len 

    ss = "wc "//trim(adjustl(filename))//" > file_length_tmp"
    call system(ss)
    open(file='file_length_tmp',action='read',unit=111)
    read(111,*) file_len
    close(111)
    return 
  end function file_len
