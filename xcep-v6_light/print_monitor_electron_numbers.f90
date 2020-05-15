
  subroutine print_monitor_electron_numbers(file_dfet_out,nspin,nelectr)
     implicit none 
     integer :: nspin, file_dfet_out
     real(8) :: nelectr(nspin)
     !
     if (nspin==1) write(file_dfet_out,'(a,f12.6)')'total electron number: ',nelectr(1)
     if (nspin==2) then 
        write(file_dfet_out,'(a,f12.6)')'total spin-up electron number:   ',nelectr(1)
        write(file_dfet_out,'(a,f12.6)')'total spin-down electron number: ',nelectr(2)
        write(file_dfet_out,'(a,f12.6)')'total electron: ',sum(nelectr)
        write(file_dfet_out,'(a,f12.6)')'mag:            ',nelectr(1)-nelectr(2)
     endif
  end subroutine 
