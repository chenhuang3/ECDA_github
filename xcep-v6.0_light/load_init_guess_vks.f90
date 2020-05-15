!
! load initial KS potential
!
subroutine load_init_guess_vks(nfft,nspin,file_dfet_out,global_vks)

  implicit none 
  integer :: nfft,nspin,isp,ii,file_dfet_out
  real(8) :: global_vks(nfft,nspin)
  !
  open(file='init_guess_total_vks.dat',unit=111,action='read')
  do isp=1,nspin
    do ii=1,nfft
      read(111,*)global_vks(ii,isp)
    enddo
  enddo
  close(111)

end subroutine load_init_guess_vks


