!
! remove the means of each spin channel
!
subroutine remove_mean(nspin,nfft,arr)

  implicit none 
  integer :: nspin,nfft,isp
  real(8) :: arr(nfft,nspin) 

  do isp=1,nspin 
    arr(:,isp) = arr(:,isp) - sum(arr(:,isp))/nfft
  enddo

endsubroutine remove_mean 
