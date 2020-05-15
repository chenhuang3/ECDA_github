
!==========================================
! Compute the gradient of arr(nfft)
!
! Chen Huang Feb/2011
!==========================================
subroutine gradient(nfft,n1,n2,n3,arr,qvec,g)
  implicit none

  integer, intent(in) :: nfft,n1,n2,n3
  real(kind=8),intent(in)   :: arr(nfft),qvec(3,(n1/2+1)*n2*n3)
  real(kind=8),intent(out)  :: g(3,nfft)
 
  ! local vars
  !===========
  integer         :: dim1,kk
  real(kind=8)    :: q3d(3,n1/2+1,n2,n3),rtmp(n1,n2,n3)
  complex(kind=8) :: fft1((n1/2+1),n2,n3),fft2((n1/2+1),n2,n3)


  !! function begins !!
  !=====================
  
  call fft(n1,n2,n3,reshape(arr,(/n1,n2,n3/)),fft1,1)

  dim1 = (n1/2+1)

  q3d(1,:,:,:) = reshape(qvec(1,:),(/dim1,n2,n3/))
  q3d(2,:,:,:) = reshape(qvec(2,:),(/dim1,n2,n3/))
  q3d(3,:,:,:) = reshape(qvec(3,:),(/dim1,n2,n3/))

  do kk=1,3 
    fft2 = fft1*dcmplx(0.d0,q3d(kk,:,:,:))
    call fft(n1,n2,n3,rtmp,fft2,-1)
    g(kk,:) = reshape(rtmp,(/nfft/))
  enddo

  return
end subroutine gradient
