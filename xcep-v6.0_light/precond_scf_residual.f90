!
! The code follows the idea in paper:
!
! "Preconditioning of self-consistent-field cycles in density-functional theory: The extrapolar method"
! Anglade and Gonze (PRB 78, 045126 2008)
! 
! See equations (14) and (15) in above paper.
! The model dielectric function follows ABINIT 
! https://docs.abinit.org/variables/gstate/#dielng
! 
!                 1 + L^2 * q^2
!  diel(q^2)= ------------------------
!              (1/diemac + L^2 * q^2)
!
!  Note that: at q=0, diel(q) is imposed to 1. 
!
subroutine precond_scf_residual(n1,n2,n3,nfft,qvec,diemac,dielen,x,Px)
 implicit none 

 integer :: n1,n2,n3,nfft, dim1, ix, iy, iz
 real(8) :: x(nfft), Px(nfft), & 
            q3d(3,n1/2+1,n2,n3), & 
            qvec(3,(n1/2+1)*n2*n3), & 
            factor_k, q0, & 
            rtmp(n1,n2,n3), & 
            qq,die,diemac,dielen

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

       ! make model dielectric function 
       if (qq>1e-12) then 
         die = (1.d0+dielen**2*qq)/(1.d0/diemac+dielen**2*qq)
       else
         die = 1.d0 ! impose to 1 at q=0 (the potential will not shifted in space)
       endif 

       fft2(ix,iy,iz) = fft1(ix,iy,iz)/die
     enddo
   enddo
 enddo

 call FFT(n1,n2,n3,rtmp,fft2,-1)
 Px = reshape(rtmp,(/n1*n2*n3/))
end subroutine 

