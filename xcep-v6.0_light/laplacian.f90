!==========================================
! Compute the Laplacian of arr(nfft)
!
! Chen Huang Feb/2011
!==========================================
subroutine laplacian(nfft,n1,n2,n3,arr,qvec,lap)
  implicit none

  integer, intent(in) :: nfft,n1,n2,n3
  real(kind=8),intent(in)   :: arr(nfft),qvec(3,(n1/2+1)*n2*n3)
  real(kind=8),intent(out)  :: lap(nfft)
 
  ! local vars
  !===========
  integer         :: i,dim1,index,kk
  complex(kind=8) :: fft1((n1/2+1),n2,n3)
  real(kind=8)    :: qnorm(n1/2+1,n2,n3),tmp(n1,n2,n3)
 
  call fft(n1,n2,n3,reshape(arr,(/n1,n2,n3/)),fft1,1)

  qnorm = reshape(sum(qvec**2,1),(/n1/2+1,n2,n3/))

  fft1 = -fft1*qnorm

  call fft(n1,n2,n3,tmp,fft1,-1)

  lap = reshape(tmp,(/nfft/))

  return
end subroutine laplacian
