 !===========================================================
 ! precondition x
 !
 subroutine precond_kerker(n1,n2,n3,nfft,qvec,q0,x,Px)
  implicit none 

  integer :: n1,n2,n3,nfft, dim1, ix, iy, iz
  real(8) :: x(nfft), Px(nfft), & 
             q3d(3,n1/2+1,n2,n3), & 
             qvec(3,(n1/2+1)*n2*n3), & 
             factor_k, q0, & 
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

        fft2(ix,iy,iz) = fft1(ix,iy,iz)*qq/(qq+q0)
      enddo
    enddo
  enddo

  call FFT(n1,n2,n3,rtmp,fft2,-1)
  Px = reshape(rtmp,(/n1*n2*n3/))
 end subroutine 






 !===========================================================
 ! precondition x
 !
 subroutine precond_kerker_dfet(n1,n2,n3,nfft,qvec,q0,x,Px)
  implicit none 

  integer :: n1,n2,n3,nfft, dim1, ix, iy, iz
  real(8) :: x(nfft), Px(nfft), & 
             q3d(3,n1/2+1,n2,n3), & 
             qvec(3,(n1/2+1)*n2*n3), & 
             factor_k, q0,  qq_min, & 
             rtmp(n1,n2,n3), metric, & 
             qq((n1/2+1),n2,n3),qq1

  complex(kind=8) :: fft1(n1/2+1,n2,n3), & 
                     fft2(n1/2+1,n2,n3)

  dim1=(n1/2+1)
  call FFT(n1,n2,n3,reshape(x,(/n1,n2,n3/)),fft1,1)

  ! convert qvec to 3-dimension 
  q3d(1,:,:,:) = reshape(qvec(1,:),(/dim1,n2,n3/))
  q3d(2,:,:,:) = reshape(qvec(2,:),(/dim1,n2,n3/))
  q3d(3,:,:,:) = reshape(qvec(3,:),(/dim1,n2,n3/))

  qq = sum(q3d**2,1)

  ! get the smallest nonzero qq
  qq_min = 1e10
  do ix=1,dim1
    do iy=1,n2
      do iz=1,n3
        qq1 = qq(ix,iy,iz)
        if ( qq1>1e-8 .and. qq_min>qq1 ) then 
          qq_min = qq1
        endif 
      enddo
    enddo
  enddo


  ! apply inverse of kerker precond
  do ix=1,dim1
    do iy=1,n2
      do iz=1,n3
        qq1 = qq(ix,iy,iz)

        if ( qq1>1e-8 ) then 
          metric = (qq1+q0)/qq1
        else 
          metric = (qq_min+q0)/qq_min
        endif 
        fft2(ix,iy,iz) = fft1(ix,iy,iz) * 1.d0/metric**2
      enddo
    enddo
  enddo

  call FFT(n1,n2,n3,rtmp,fft2,-1)
  Px = reshape(rtmp,(/n1*n2*n3/))
 end subroutine 




 !===========================================================
 ! precondition x
 !
 subroutine precond_kerker_inv(n1,n2,n3,nfft,qvec,q0,x,Px)
  implicit none 

  integer :: n1,n2,n3,nfft, dim1, ix, iy, iz
  real(8) :: x(nfft), Px(nfft), & 
             q3d(3,n1/2+1,n2,n3), & 
             qvec(3,(n1/2+1)*n2*n3), & 
             factor_k, q0,  qq_min, & 
             rtmp(n1,n2,n3), & 
             qq((n1/2+1),n2,n3),qq1

  complex(kind=8) :: fft1(n1/2+1,n2,n3), & 
                     fft2(n1/2+1,n2,n3)

  dim1=(n1/2+1)
  call FFT(n1,n2,n3,reshape(x,(/n1,n2,n3/)),fft1,1)

  ! convert qvec to 3-dimension 
  q3d(1,:,:,:) = reshape(qvec(1,:),(/dim1,n2,n3/))
  q3d(2,:,:,:) = reshape(qvec(2,:),(/dim1,n2,n3/))
  q3d(3,:,:,:) = reshape(qvec(3,:),(/dim1,n2,n3/))

  qq = sum(q3d**2,1)

  ! get the smallest nonzero qq
  qq_min = 1e10
  do ix=1,dim1
    do iy=1,n2
      do iz=1,n3
        qq1 = qq(ix,iy,iz)
        if ( qq1>1e-8 .and. qq_min>qq1 ) then 
          qq_min = qq1
        endif 
      enddo
    enddo
  enddo

  ! apply inverse of kerker precond
  do ix=1,dim1
    do iy=1,n2
      do iz=1,n3
        qq1 = qq(ix,iy,iz)
        if ( qq1>1e-8 ) then 
           fft2(ix,iy,iz) = fft1(ix,iy,iz) * (qq1+q0)/qq1
        else 
           fft2(ix,iy,iz) = fft1(ix,iy,iz) * (qq_min+q0)/qq_min
        endif 
      enddo
    enddo
  enddo

  call FFT(n1,n2,n3,rtmp,fft2,-1)
  Px = reshape(rtmp,(/n1*n2*n3/))

 end subroutine  precond_kerker_inv

