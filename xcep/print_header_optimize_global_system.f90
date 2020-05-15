  subroutine print_header_optimize_global_system(file_dfet_out,iter_scf) 
   implicit none 
   integer :: iter_scf,file_dfet_out
   print *,''
   print *,''
   print *,'                    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
   write(*,'(a,i3)')"   optimizing global vks, iter: ",iter_scf
   print *,'                    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
   print *,''
   print *,''
   write(file_dfet_out,'(a)')''
   write(file_dfet_out,'(a)')'   ===================================='
   write(file_dfet_out,'(a,i3)')"   optimizing global vks, iter: ",iter_scf
   write(file_dfet_out,'(a)')'   ===================================='
   write(file_dfet_out,'(a)')''
  end subroutine print_header_optimize_global_system
