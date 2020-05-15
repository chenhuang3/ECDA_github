!
! compute Hess|x>, where 
!
! Hess = d^2E/drho(r)drho(r')
!
! Created by Chen Huang 12/27/2019
!
subroutine apply_hess_vwreg1(nfft,n1,n2,n3,qvec,rho,x,hx)

 use comm_data 

 implicit none 

 integer :: nfft, n1, n2, n3

 real(8),parameter :: pi=3.1415926d0, & 
                      cTF = 3.d0/10.d0*(3.d0*pi**2)**(2.d0/3.d0)

 real(8) :: rho(nfft), & 
            qvec(3,(n1/2+1)*n2*n3), & 
            rho_reg(nfft), &
            lap1(nfft), & 
            lap2(nfft), & 
            x(nfft), hx(nfft), & 
            tf_kernel(nfft)  ! TF kernel 


 !>>>>>>>>>> function begins <<<<<<<<<<<<<<

 rho_reg = rho + vw_reg
 
 tf_kernel = 5.d0/3.d0*cTF*2.d0/3.d0*rho**(-1.d0/3.d0) 

 call laplacian(nfft,n1,n2,n3,  sqrt(rho_reg),qvec,lap1)
 call laplacian(nfft,n1,n2,n3,x/sqrt(rho_reg),qvec,lap2)

 hx =   tf_kernel*x & 
      + vw_lam*0.25d0*lap1/rho_reg**1.5d0*x & 
      - vw_lam*0.25d0/sqrt(rho_reg)*lap2

end subroutine
