 !
 ! precondition x
 !
 subroutine precond_teter_inv(n1,n2,n3,nfft,qvec,ke,x,Px)
  implicit none 

  integer :: n1,n2,n3,nfft, dim1, ix, iy, iz
  real(8) :: x(nfft), Px(nfft), & 
             q3d(3,n1/2+1,n2,n3), & 
             qvec(3,(n1/2+1)*n2*n3), & 
             factor_k, ke, & 
             rtmp(n1,n2,n3), & 
             qq

  complex(kind=8) :: fft1(n1/2+1,n2,n3), & 
                     fft2(n1/2+1,n2,n3)

  dim1=(n1/2+1)
  call FFT(n1,n2,n3,reshape(x,(/n1,n2,n3/)),fft1,1)

  ! convert qvec to 3-dimension 
  q3d(1,:,:,:) = reshape(qvec(1,:),(/dim1,n2,n3/))
  q3d(2,:,:,:) = reshape(qvec(2,:),(/dim1,n2,n3/))
  q3d(3,:,:,:) = reshape(qvec(3,:),(/dim1,n2,n3/))


  ! make the 4*pi/q^2*rhor(q) and store in fft2(q)
  do ix=1,dim1
    do iy=1,n2
      do iz=1,n3
        qq = q3d(1,ix,iy,iz)**2 + & 
             q3d(2,ix,iy,iz)**2 + & 
             q3d(3,ix,iy,iz)**2

        qq = qq/2.d0/ke ! x in Payne's paper in Eq.5.16
        factor_k = (27.d0+18.d0*qq+12.d0*qq*qq+8.d0*qq**3)/ & 
                   (27.d0+18.d0*qq+12.d0*qq*qq+8.d0*qq**3+16.d0*qq**4)

        fft2(ix,iy,iz) = fft1(ix,iy,iz)/factor_K
      enddo
    enddo
  enddo

  call FFT(n1,n2,n3,rtmp,fft2,-1)
  Px = reshape(rtmp,(/n1*n2*n3/))

 end subroutine 
