
subroutine print_optimize_w_header(file_dfet_out,nfft,nspin,gmax_limit,cgw_tol,u)

  integer :: file_dfet_out,nfft,nspin
  real(8) :: gmax_limit, cgw_tol, u(nfft,nspin)

  if (nspin==2) then 
    write(file_dfet_out,'(a,2es12.4)')'(initial embpot) u_alpha: ',minval(u(:,1)),maxval(u(:,1))
    write(file_dfet_out,'(a,2es12.4)')'                 u_beta:  ',minval(u(:,2)),maxval(u(:,2))
  else 
    write(file_dfet_out,'(a,2es12.4)')'(initial embpot) u: ',minval(u(:,1)),maxval(u(:,1))
  endif 
  !
  print *,''
  print *,' -------------------------------------------------------------------- '
  print *,'   Optimizing embedding potential for fixed subsystem KS potentials  '
  print *,' -------------------------------------------------------------------- '
  print *,''
  write(file_dfet_out,'(a)')'minimizing W functional ...'
  write(file_dfet_out,'(a,es12.4)')'gmax_limit: ',gmax_limit
  write(file_dfet_out,'(a,es12.4)')'cgw_tol:    ',cgw_tol
  call flush(file_dfet_out)

end subroutine print_optimize_w_header       
