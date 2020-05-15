! 
! for given density, compute the hartree energy and potential 
! the spin density will be sumed to yield total density first 
!
! The Coulumb potential is screended by error function as in the case of HSE 
! see the paper by Kresse JCP 124, 154709 (2006) Eq. (10)
!
subroutine hartree_screened(n1,n2,n3,nspin,ucvol,qvec,rhor,vhart,ehart)

  implicit none 

  integer :: n1,n2,n3,ix,iy,iz,dim1,kk,isp,nspin

  real(kind=8) :: rhor(n1*n2*n3,nspin), &   ! spin density
                  vhart(n1*n2*n3), &        ! hartree potential due to rho_up + rho_dn
                  qvec(3,(n1/2+1)*n2*n3), &
                  rtmp(n1,n2,n3), &         ! temporary array in real space
                  ucvol,pi,     &           ! volume of cell
                  q3d(3,n1/2+1,n2,n3), mu, &
                  qq, &                     ! q^2
                  ehart                     ! hartree energy
  
  complex(kind=8) :: fft1(n1/2+1,n2,n3), & 
                     fft2(n1/2+1,n2,n3)

  dim1=(n1/2+1)

  call FFT(n1,n2,n3,reshape(sum(rhor,2),(/n1,n2,n3/)),fft1,1)

  ! convert qvec to 3-dimension 
  q3d(1,:,:,:) = reshape(qvec(1,:),(/dim1,n2,n3/))
  q3d(2,:,:,:) = reshape(qvec(2,:),(/dim1,n2,n3/))
  q3d(3,:,:,:) = reshape(qvec(3,:),(/dim1,n2,n3/))

  vhart = 0.0d0

  ! make the 1/q^2*rhor(q) and store in fft2(q)
  do ix=1,dim1
    do iy=1,n2
      do iz=1,n3
        qq = q3d(1,ix,iy,iz)**2 + & 
             q3d(2,ix,iy,iz)**2 + & 
             q3d(3,ix,iy,iz)**2
        !!qq = qq/4.d0/3.1415926d0
        if (ix==1 .and. iy==1 .and. iz==1) then
          fft2(ix,iy,iz) = (0.0D0,0.0D0)
        else
          mu = (1/3.0d0) * 0.5291772106d0 ! 1/3 angstrom^-1 -> bohr^-1
          pi = 3.1415926d0
          fft2(ix,iy,iz) = 4.d0*pi*fft1(ix,iy,iz)/qq*(1.d0-exp(-qq/4.d0/mu**2))
        endif
      enddo
    enddo
  enddo

  call FFT(n1,n2,n3,rtmp,fft2,-1)
  vhart = reshape(rtmp,(/n1*n2*n3/))

  ! get hartree energy 
  ehart = sum(vhart*sum(rhor,2))*ucvol/dble(n1*n2*n3)*0.5D0

  !print *, ''
  !print *, '------ hartree.f90 --------'
  !write(6,'(a,f12.6,a,2f12.4)') "  ehart: ",ehart, & 
  !' Ha. min/max(vhart):',minval(vhart),maxval(vhart)
  !print *, ''

end subroutine hartree_screened
