  subroutine check_file_length(nspin,n1,n2,n3)
    implicit none 
    integer :: nspin,n1,n2,n3
    integer :: file_len 
    ! check the length of the file 
    if (nspin==2) then 
   !   if ( file_len('init_guess_total_rho.dat')  /=  2*n1*n2*n3) then 
   !     print *,'init_guess_total_rho.dat length is not 2*n1*n2*n3, error stop!'
   !     stop 
   !   endif 
      if ( file_len('init_guess_total_vks.dat')  /=  2*n1*n2*n3) then 
        print *,'init_guess_total_vks.dat length is not 2*n1*n2*n3, error stop!'
        call flush(6)
        stop 
      endif
    endif 
    if (nspin==1) then 
   !   if ( file_len('init_guess_total_rho.dat')  /=  1*n1*n2*n3) then 
   !     print *,'init_guess_total_rho.dat length is not 2*n1*n2*n3, error stop!'
   !     stop 
   !   endif 
      if ( file_len('init_guess_total_vks.dat')  /=  1*n1*n2*n3) then 
        print *,'init_guess_total_vks.dat length is not 2*n1*n2*n3, error stop!'
        call flush(6)
        stop 
      endif
    endif 
  end subroutine check_file_length
