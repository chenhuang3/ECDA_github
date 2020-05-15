!
! Precondition a vector in the Fourier space with the matrix
!
!   Q = 1/(q^2/(q^2+aa) + reg_param)
!
! On exit, v is updated 
!
subroutine precond(n1,n2,n3,ucvol,reg,qvec,v)

  implicit none 

  integer :: n1,n2,n3,ix,iy,iz,dim1,kk,isp,nspin

  real(kind=8) :: v(n1*n2*n3), &   ! spin density
                  reg, & 
                  qvec(3,(n1/2+1)*n2*n3), &
                  rtmp(n1,n2,n3), &         ! temporary array in real space
                  ucvol,        &           ! volume of cell
                  q3d(3,n1/2+1,n2,n3), &
                  qq, &                     ! q^2
                  ehart                     ! hartree energy
  
  complex(kind=8) :: fft1(n1/2+1,n2,n3), & 
                     fft2(n1/2+1,n2,n3)


  ! local vars 
  real(8) :: vtmp(n1,n2,n3),diemac,dielen,chi_q,cc,pi=3.1415926d0
  
  diemac = 2.0
  dielen = 1.0774841d0


  dim1=(n1/2+1)
  call FFT(n1,n2,n3,reshape(v,(/n1,n2,n3/)),fft1,1)

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

        !!chi_q = qq/4.0d0/pi*(1.d0-1.d0/diemac)/(1.d0/diemac+qq*dielen**2)
        chi_q = qq/(qq+0.1)
        cc = reg
        !chi_q = qq/(qq+1.0)

        fft2(ix,iy,iz) = fft1(ix,iy,iz)/(chi_q+cc)
      enddo
    enddo
  enddo

  call FFT(n1,n2,n3,rtmp,fft2,-1)
  v = reshape(rtmp,(/n1*n2*n3/))

end subroutine precond
