subroutine check_subsystem_charges(dvol,nspin,nfft,nsys,ref_rho,sub_rhor)
  implicit none 
  integer :: isp, nspin, nsys, nfft
  real(8) :: ref_rho(nfft,nspin), sub_rhor(nfft,nspin,nsys)
  real(8) :: totalQ, totalQ1, dvol
  
  ! check if the total charge is equal to the sum of subsystem charges 
  if (nspin==1) then 
    totalQ = sum(ref_rho)*dvol 
    totalQ1 = sum(sub_rhor)*dvol
    if ( abs(totalQ - totalQ1)>1e-4) then 
      print *,'totalQ_ref: ',totalQ
      print *,'totalQ_subsystem: ',totalQ1,' does not match, stop'
      stop
    endif
  endif
  if (nspin==2) then 
    do isp=1,nspin
    totalQ = sum(ref_rho(:,isp))*dvol
    totalQ1 = sum(sub_rhor(:,isp,:))*dvol
    if ( abs(totalQ - totalQ1)>1e-4) then 
      print *,'totalQ_ref: ',totalQ
      print *,'totalQ_subsystem: ',totalQ1,' does not match, stop'
      stop
    endif
    enddo
  endif
end subroutine 
