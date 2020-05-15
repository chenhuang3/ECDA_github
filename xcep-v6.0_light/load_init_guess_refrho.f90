  !
  ! ref_rho might be from L(S)DA calculation of the total system
  !
  subroutine load_init_guess_refrho(nfft, nspin, dvol, file_dfet_out, ref_rho)
    implicit none 
    integer :: nfft,nspin,ii,isp,file_dfet_out
    real(8) :: ref_rho(nfft,nspin),dvol
    !
    open(file='init_guess_total_rho.dat',unit=111,action='read',form='formatted')
    do isp=1,nspin
      print *,'loading spin channel: ',isp,' ...'
      do ii=1,nfft 
        read(111,*)ref_rho(ii,isp)
      enddo
    enddo
    close(111)
    write(file_dfet_out,'(a)')'File init_guess_total_rho.dat is loaded and set to ref_rho.'
    if (nspin==2) then 
      print *,'total_Q (alpha): ',sum(ref_rho(:,1))*dvol
      print *,'total_Q (beta ): ',sum(ref_rho(:,2))*dvol
    else 
      print *,'total_Q: ',sum(ref_rho(:,1))*dvol
    endif
    print *,''
  end subroutine load_init_guess_refrho

