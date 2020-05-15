subroutine get_reg_matrix(nspin,nfft,rho,reg_mat)

 implicit none 
 integer :: nfft,nspin 
 real(8) :: eps, rho0, sigma, reg_mat(nfft), rho(nfft,nspin)

 eps = 1.0;
 rho0 = 1e-4;
 sigma = rho0*1.0; 

! y = eps*(1+exp(-rho0/sigma))./(1+exp((xx-rho0)/sigma));
! figure
! semilogx(xx,y)

 reg_mat = eps*(1.d0+exp(-rho0/sigma))/(1.d0+exp((sum(rho,2)-rho0)/sigma));

end subroutine 

