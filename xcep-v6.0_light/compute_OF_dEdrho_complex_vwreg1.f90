subroutine compute_OF_dEdrho_complex_vwreg1(nfft,n1,n2,n3,ucvol,qvec,rho,dEdrho)
   use comm_data 

   implicit none 
   integer    :: nfft, n1,n2,n3
   complex(8) :: rho(nfft),  & 
                 qvec(3,(n1/2+1)*n2*n3), & 
                 dEdrho(nfft),  &
                 sq_rhor(nfft), & 
                 tfpot(nfft), vext(nfft),  & 
                 lap(nfft), vw_pot(nfft)

   real(kind=8) :: ucvol, & 
                   ctf = 2.87123400018819d0, & 
                   lap_re(nfft), & 
                   lap_im(nfft)

   sq_rhor = sqrt(rho + vw_reg)   ! regularized 

   call laplacian(nfft,n1,n2,n3, real(sq_rhor),qvec,lap_re)
   call laplacian(nfft,n1,n2,n3,aimag(sq_rhor),qvec,lap_im)
   lap = cmplx(lap_re,lap_im)

   vw_pot = - 0.5d0*lap/sq_rhor
   tfpot = ctf * rho**(2.0d0/3.d0) * (5.d0/3.d0)
   dEdrho = tfpot + vw_pot*vw_lam + vext
end subroutine 

